// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2020, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

//
// ustacks -- build denovo stacks
//

#include "ustacks.h"
#include "log_utils.h"

using namespace std;

//
// Global variables to hold command-line options.
//
FileT   in_file_type;
string  in_file;
string  prefix_path;
int     num_threads       = 1;
int     sample_id         = -1;
bool    call_sec_hapl     = true;
bool    set_kmer_len      = true;
int     global_kmer_len   = 0;
int     min_merge_cov     = 3;
uint    max_subgraph      = 3;
int     dump_graph        = 0;
int     retain_rem_reads  = false;
int     deleverage_stacks = 0;
bool    remove_rep_stacks = true;
double  removal_threshold = 3.0;
int     max_utag_dist     = 2;
int     max_rem_dist      = -1;
bool    force_diff_len    = false;
bool    gapped_alignments = true;
double  min_match_len     = 0.80;
double  max_gaps          = 2.0;
int     removal_trigger;

//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.05;

int main (int argc, char* argv[]) {
try{
    parse_command_line(argc, argv);

    //
    // Set the max remainder distance to be greater than the max_utag_dist, if it is not
    // specified on the command line.
    //
    if (max_rem_dist == -1) max_rem_dist = max_utag_dist + 2;

    cerr << "ustacks parameters selected:\n"
         << "  Input file: '" << in_file   << "'\n"
         << "  Sample ID: "   << sample_id << "\n"
         << "  Min depth of coverage to create a stack (m): " << min_merge_cov << "\n"
         << "  Repeat removal algorithm: " << (remove_rep_stacks ? "enabled" : "disabled") << "\n"
         << "  Max distance allowed between stacks (M): " << max_utag_dist << "\n"
         << "  Max distance allowed to align secondary reads: " << max_rem_dist << "\n"
         << "  Max number of stacks allowed per de novo locus: " << max_subgraph << "\n"
         << "  Deleveraging algorithm: " << (deleverage_stacks ? "enabled" : "disabled") << "\n";
    if (gapped_alignments) {
        cerr << "  Gapped assembly: enabled\n"
             << "  Minimum alignment length: " << min_match_len << "\n";
    } else {
        cerr << "  Gapped assembly: disabled\n";
    }
    cerr << "  Model type: ";
    switch (model_type) {
    case snp:
        cerr << "SNP\n";
        break;
    case ::fixed:
        cerr << "Fixed\n";
        break;
    case bounded:
        cerr << "Bounded; lower epsilon bound: " << bound_low << "; upper bound: " << bound_high << "\n";
        break;
    default:
        DOES_NOT_HAPPEN;
        break;
    }
    cerr << "  Alpha significance level for model: " << alpha << "\n";
    if (force_diff_len)
        cerr << "  Forcing the allowance of sequences of different length.\n";
    cerr << flush;

    //
    // Set limits to call het or homozygote according to chi-square distribution with one
    // degree of freedom:
    //   http://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_.CF.872_value_vs_p-value
    //
    set_model_thresholds(alpha);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Load the reads and build primary & secondary stacks.
    //

    cerr << "\n" << "Loading RAD-Tags...";

    size_t n_reads, n_u_reads, n_r_reads, n_rr_reads;
    DNASeqHashMap radtags;
    map<int, Stack *> unique;
    map<int, Rem *>   remainders;
    load_radtags(in_file, radtags, n_reads);
    reduce_radtags(radtags, unique, remainders, n_u_reads, n_r_reads);
    radtags = DNASeqHashMap();

    cerr << "\n"
         << "Loaded " << n_reads << " reads; formed:\n"
         << "  " << unique.size() << " stacks representing "
         << n_u_reads << " primary reads (" << as_percentage((double)n_u_reads/n_reads) << ")\n"
         << "  " << remainders.size() << " secondary stacks representing "
         << n_r_reads << " secondary reads (" << as_percentage((double)n_r_reads/n_reads) << ")\n"
         << flush;

    //
    // Initialize the MergedStack object.
    //
    map<int, MergedStack *> merged;
    populate_merged_tags(unique, merged);

    double cov_mean, cov_stdev, cov_max, n_used_reads;
    ostream cerr_copy (cerr.rdbuf());
    cerr_copy << std::fixed << std::setprecision(2);
    auto report_cov = [&]() {
        cerr_copy << "mean=" << cov_mean << "; stdev=" << cov_stdev << "; max=" << size_t(cov_max)
                  << "; n_reads=" << size_t(n_used_reads) << "(" << as_percentage(n_used_reads/n_reads) << ")";
    };

    calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
    cerr << "\n";
    cerr << "Stack coverage: ";
    report_cov();
    cerr << "\n";

    //
    // Remove highly repetitive stacks.
    //
    size_t n_high_cov = 0;
    if (remove_rep_stacks) {
        removal_trigger = (int) floor(cov_mean + cov_stdev * removal_threshold + 1);
        cerr << "Removing repetitive stacks: cov > " << removal_trigger << " (mean+" << removal_threshold << "*stdev)...\n";
        calc_kmer_distance(merged, 1);
        assert(merged.size() == unique.size());
        n_high_cov = remove_repetitive_stacks(merged);
        cerr << "  Blacklisted " << n_high_cov << " stacks.\n";

        calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
        cerr << "Coverage after repeat removal: ";
        report_cov();
        cerr << "\n\n";
    }

    //
    // Assemble loci (merge primary stacks).
    //
    cerr << "Assembling stacks (max. dist. M=" << max_utag_dist << ")...\n";
    calc_kmer_distance(merged, max_utag_dist);
    size_t n_blacklisted;
    merge_stacks(merged, n_blacklisted);
    call_consensus(merged, unique, remainders, false);
    cerr << "  Assembled " << unique.size()-n_high_cov << " stacks into " << merged.size()
         << "; blacklisted " << n_blacklisted << " stacks.\n";
    calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
    cerr << "Coverage after assembling stacks: ";
    report_cov();
    cerr << "\n\n";

    //
    // Merge secondary stacks.
    //
    cerr << "Merging secondary stacks (max. dist. N=" << max_rem_dist << " from consensus)...\n";
    size_t utilized    = merge_remainders(merged, unique, remainders);
    size_t util_gapped = merge_gapped_remainders(merged, unique, remainders);
    call_consensus(merged, unique, remainders, false);
    cerr << "  Merged " << utilized+util_gapped << " out of " << n_r_reads << " secondary reads ("
         << as_percentage( (double) (utilized+util_gapped) / (double)n_r_reads) << "), "
         << util_gapped << " merged with gapped alignments.\n";
    n_rr_reads = n_r_reads - utilized - util_gapped;

    calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
    cerr << "Coverage after merging secondary stacks: ";
    report_cov();
    cerr << "\n\n";

    //
    // Merge loci based on alignments.
    //
    if (gapped_alignments) {
        cerr << "Assembling stacks, allowing for gaps (min. match length " << as_percentage(min_match_len) << ")...\n";
        const size_t n_ungapped_loci = merged.size();
        search_for_gaps(merged);
        merge_gapped_alns(unique, remainders, merged);
        call_consensus(merged, unique, remainders, false);
        cerr << "  Assembled " << n_ungapped_loci << " stacks into " << merged.size() << " stacks.\n";

        calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
        cerr << "Coverage after gapped assembly: ";
        report_cov();
        cerr << "\n\n";

        cerr << "Merging secondary stacks, allowing for gaps (min. match length " << as_percentage(min_match_len) << ")...\n";
        search_for_gapped_remainders(merged, unique, remainders);
        util_gapped = merge_gapped_remainders(merged, unique, remainders);
        cerr << "  Merged " << util_gapped << " out of " << n_rr_reads << " secondary reads ("
             << as_percentage( (double)util_gapped / (double)n_rr_reads) << ").\n";

        calc_coverage_distribution(merged, cov_mean, cov_stdev, cov_max, n_used_reads);
        cerr << "Coverage after merging gapped secondary stacks: ";
        report_cov();
        cerr << "\n\n";
    }

    //
    // Report how many reads were used.
    //
    cerr << "Final coverage: ";
    report_cov();
    cerr << "\n";

    //
    // Call the final consensus sequence and invoke the SNP model.
    //
    cerr << "Calling consensus sequences and haplotypes for catalog assembly...\n";
    call_consensus(merged, unique, remainders, true);

    //
    // Write output files.
    //
    cerr << "Writing tags, SNPs, and alleles files...\n";
    write_results(merged, unique, remainders);
    cerr << "done.\n";

    cerr << "ustacks is done.\n";
    return 0;
} catch(std::exception& e) {
    return stacks_handle_exceptions(e);
}}

int
merge_gapped_alns(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged)
{
    map<int, MergedStack *> new_merged;
    MergedStack *tag_1, *tag_2, *merged_tag;

    int  id        = 1;
    uint merge_cnt = 0;

    set<int> processed;
    string   cigar_1, cigar_2;
    vector<pair<char, uint> > cigar;

    for (auto it = merged.begin(); it != merged.end(); it++) {
        if (processed.count(it->first) > 0)
            continue;

        tag_1 = it->second;
        sort(tag_1->alns.begin(), tag_1->alns.end(), rank_alignments);

        //
        // No gapped alignments, or no optimal alignment for this stack, or
        // this stack has already been set aside.
        //
        if (tag_1->masked || tag_1->alns.size() == 0)
            continue;
        if (tag_1->alns.size() > 1 && tag_1->alns[0].pct_id == tag_1->alns[1].pct_id)
            continue;

        //
        // Found one or more gapped alignments. Make sure the best alignment for each of the aligned pairs
        // of tags are reciprocal to each other.
        //
        tag_2 = merged[tag_1->alns[0].id];
        sort(tag_2->alns.begin(), tag_2->alns.end(), rank_alignments);

        if (tag_2->masked || tag_2->alns.size() == 0)
            continue;
        if (tag_2->alns.size() > 1 && tag_2->alns[0].pct_id == tag_2->alns[1].pct_id)
            continue;

        if (tag_1->id != (int)tag_2->alns[0].id)
            continue;

        cigar_1 = invert_cigar(tag_1->alns[0].cigar);
        cigar_2 = tag_2->alns[0].cigar;

        if (cigar_1 == cigar_2) {
            parse_cigar(tag_1->alns[0].cigar.c_str(), cigar);

            //
            // Check that the alignment still contains fewer than
            // max_utag_dist mismatches.
            //
            if (dist(tag_1->con, tag_2->con, cigar) > max_utag_dist)
                continue;

            //
            // If the alignment has too many gaps, skip it.
            //
            if (tag_1->alns[0].gap_cnt > (max_gaps + 1))
                continue;

            //
            // If the alignment doesn't span enough of the two sequences, skip it.
            //
            if (tag_1->alns[0].pct_id < min_match_len)
                continue;

            //
            // Edit the sequences to accommodate any added gaps.
            //
            edit_gapped_seqs(unique, rem, tag_1, cigar);

            parse_cigar(tag_2->alns[0].cigar.c_str(), cigar);
            edit_gapped_seqs(unique, rem, tag_2, cigar);

            //
            // Merge the tags.
            //
            merged_tag     = merge_tags(tag_1, tag_2, id);
            new_merged[id] = merged_tag;
            id++;

            //
            // Record the gaps.
            //
            uint pos = 0;
            for (uint j = 0; j < cigar.size(); j++) {
                if (cigar[j].first == 'I' || cigar[j].first == 'D')
                    merged_tag->gaps.push_back(Gap(pos, pos + cigar[j].second));
                pos += cigar[j].second;
            }

            processed.insert(tag_1->id);
            processed.insert(tag_2->id);

            merge_cnt++;
        }
    }

    set<int> merge_set;
    for (auto it = merged.begin(); it != merged.end(); it++) {
        if (processed.count(it->first))
            continue;
        tag_1          = it->second;
        merge_set.insert(tag_1->id);
        tag_2          = merge_tags(merged, merge_set, id);
        new_merged[id] = tag_2;
        merge_set.clear();
        id++;
    }

    //
    // Free the memory from the old map of merged tags.
    //
    for (auto it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    merged = new_merged;
    return 0;
}

bool
rank_alignments(Aln a, Aln b)
{
    return a.pct_id > b.pct_id;
}

int
edit_gapped_seqs(map<int, Stack *> &unique, map<int, Rem *> &rem, MergedStack *tag, vector<pair<char, uint> > &cigar)
{
    int    stack_id;
    Stack *s;
    Rem   *r;
    char  *buf = new char[cigar_length_padded(cigar) + 1];
    string seq;

    for (uint i = 0; i < tag->utags.size(); i++) {
        stack_id = tag->utags[i];
        s = unique[stack_id];

        s->seq->seq(buf);
        seq = apply_cigar_to_seq(buf, cigar);

        delete s->seq;
        s->seq = new DNANSeq(seq.length(), seq.c_str());
    }

    for (uint i = 0; i < tag->remtags.size(); i++) {
        stack_id = tag->remtags[i];
        r = rem[stack_id];

        r->seq->seq(buf);
        seq = apply_cigar_to_seq(buf, cigar);

        delete r->seq;
        r->seq = new DNANSeq(seq.length(), seq.c_str());
    }

    delete [] buf;

    return 0;
}

int
edit_gaps(vector<pair<char, uint> > &cigar, char *seq)
{
    char *buf;
    uint  size = cigar.size();
    char  op;
    uint  dist, bp, len, buf_len, buf_size, j, k, stop;

    len = strlen(seq);
    bp  = 0;

    buf      = new char[len + 1];
    buf_size = len + 1;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }
            break;
        case 'D':
            //
            // A deletion has occured in the read relative to the reference genome.
            // Pad the read with sufficent Ns to match the deletion, shifting the existing
            // sequence down. Trim the final length to keep the read length consistent.
            //
            k = bp >= len ? len : bp;

            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len         = strlen(buf);

            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }
            break;
        case 'I':
        case 'M':
            bp += dist;
            break;
        default:
            break;
        }
    }

    delete [] buf;

    return 0;
}

int
search_for_gaps(map<int, MergedStack *> &merged)
{
    //
    // Search for loci that can be merged with a gapped alignment.
    //
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    MergedStack   *tag_1, *tag_2;

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (auto it = merged.begin(); it != merged.end(); it++)
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate.
    //
    size_t con_len = 0;
    for (auto miter = merged.begin(); miter != merged.end(); miter++)
        con_len = miter->second->len > con_len ? miter->second->len : con_len;
    size_t kmer_len  = set_kmer_len ? 19 : global_kmer_len;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    uint min_hits = (round((double) con_len * min_match_len) - (kmer_len * max_gaps)) - kmer_len + 1;

    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);

    #pragma omp parallel private(tag_1, tag_2)
    {
        KmerHashMap::iterator h;
        vector<char *> query_kmers;
        set<string>    uniq_kmers;
        GappedAln     *aln = new GappedAln(con_len);
        AlignRes       a;

        size_t num_kmers = con_len - kmer_len + 1;
        initialize_kmers(kmer_len, num_kmers, query_kmers);

        #pragma omp for schedule(dynamic)
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = merged[keys[i]];

            //
            // Don't compute distances for masked tags.
            //
            if (tag_1->masked) continue;

            //
            // Don't compare tags that are already at or above max_locus_stacks.
            //
            if (tag_1->utags.size() >= max_subgraph)
                continue;

            if (tag_1->len < kmer_len) continue;
            num_kmers = tag_1->len - kmer_len + 1;
            generate_kmers_lazily(tag_1->con, kmer_len, num_kmers, query_kmers);

            assert(num_kmers > 0);

            uniq_kmers.clear();
            for (uint j = 0; j < num_kmers; j++)
                uniq_kmers.insert(query_kmers[j]);

            vector<int> hits;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (auto j = uniq_kmers.begin(); j != uniq_kmers.end(); j++) {
                h = kmer_map.find(j->c_str());

                if (h != kmer_map.end())
                    for (uint k = 0; k <  h->second.size(); k++)
                        hits.push_back(h->second[k]);
            }

            //
            // Sort the vector of indexes; provides the number of hits to each allele/locus
            // and orders them largest to smallest.
            //
            sort(hits.begin(), hits.end());

            //
            // Iterate through the list of hits and collapse them down by number of kmer hits per allele.
            //
            uint hit_cnt, index, prev_id, allele_id, hits_size, stop;
            vector<pair<uint, uint>> ordered_hits;

            hits_size = hits.size();

            if (hits_size == 0)
                continue;

            prev_id = hits[0];
            index   = 0;

            do {
                hit_cnt   = 0;
                allele_id = prev_id;

                while (index < hits_size && (uint) hits[index] == prev_id) {
                    hit_cnt++;
                    index++;
                }

                if (index < hits_size)
                    prev_id = hits[index];

                // Don't compare tag_1 against itself.
                if (tag_1->id == (int) allele_id) continue;

                ordered_hits.push_back(make_pair(allele_id, hit_cnt));

            } while (index < hits_size);

            if (ordered_hits.size() == 0)
                continue;

            //
            // Process the hits from most kmer hits to least kmer hits.
            //
            sort(ordered_hits.begin(), ordered_hits.end(), compare_pair_intint);

            //
            // Align at least one locus, regardless of kmer count, or more than one locus if
            // there are sufficient k-mer hits supporting it.
            //
            stop = 1;
            for (uint j = 1; j < ordered_hits.size(); j++)
                if (ordered_hits[j].second < min_hits) {
                    stop = j;
                    break;
                }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            for (uint j = 0; j < stop; j++) {
                tag_2 = merged.at(ordered_hits[j].first);

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                //
                // Don't compare tags that are already at or above max_locus_stacks.
                //
                if (tag_2->utags.size() >= max_subgraph)
                    continue;

                aln->init(tag_1->len, tag_2->len);

                if (aln->align(tag_1->con, tag_2->con)) {
                    a = aln->result();
                    tag_1->alns.push_back(Aln(tag_2->id, a.cigar, a.pct_id, a.gap_cnt));
                }
            }
        }

        //
        // Free the k-mers we generated for this query.
        //
        for (uint j = 0; j < query_kmers.size(); j++)
            delete [] query_kmers[j];
        query_kmers.clear();

        delete aln;
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

size_t
merge_remainders(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem)
{
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    size_t max_rem_len = 0;
    for (auto it = rem.begin(); it != rem.end(); it++) {
        keys.push_back(it->first);
        max_rem_len = it->second->seq->size() > max_rem_len ? it->second->seq->size() : max_rem_len;
    }

    //
    // Calculate the number of k-mers we will generate based on the longest sequence present.
    // If kmer_len == 0, determine the optimal length for k-mers.
    //
    size_t con_len  = 0;
    for (auto miter = merged.begin(); miter != merged.end(); miter++)
        con_len = miter->second->len > con_len ? miter->second->len : con_len;
    size_t kmer_len = set_kmer_len ? determine_kmer_length(con_len, max_rem_dist) : global_kmer_len;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, max_rem_dist, con_len, set_kmer_len ? true : false);

    //cerr << "  Distance allowed between stacks: " << max_rem_dist
    //     << "; searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); "
    //     << min_hits << " k-mer hits required.\n";

    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
    size_t utilized = 0;

    #pragma omp parallel reduction(+: utilized)
    {
        KmerHashMap::iterator h;
        vector<char *> rem_kmers;
        Cigar      cigar;
        string     seq;
        char      *buf = new char[max_rem_len + 1];
        size_t     num_kmers = 0;

        #pragma omp for schedule(dynamic)
        for (uint j = 0; j < keys.size(); j++) {
            auto it = rem.find(keys[j]);
            Rem  *r = it->second;

            //
            // Generate the k-mers for this remainder sequence
            //
            r->seq->seq(buf);
            if (r->seq->size() < kmer_len) continue;
            num_kmers = r->seq->size() - kmer_len + 1;
            generate_kmers_lazily(buf, kmer_len, num_kmers, rem_kmers);

            assert(num_kmers > 0);

            map<int, int> hits;
            //
            // Lookup the occurances of each remainder k-mer in the MergedStack k-mer map
            //
            for (uint k = 0; k < num_kmers; k++) {
                h = kmer_map.find(rem_kmers[k]);

                if (h != kmer_map.end())
                    for (uint n = 0; n < h->second.size(); n++)
                        hits[h->second[n]]++;
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int> dists;
            for (auto hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {

                if (hit_it->second < min_hits)
                    continue;

                MergedStack *tag_1 = merged[hit_it->first];

                int d = dist(tag_1, buf);

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance
                //
                if (d <= max_rem_dist) {
                    dists[hit_it->first] = d;
                }
            }

            //
            // Check to see if there is a uniquely low distance, if so,
            // merge this remainder tag. If not, discard it, since we
            // can't locate a single best-fitting Stack to merge it into.
            //
            map<int, int>::iterator s;
            int min_id = -1;
            int count  =  0;
            int dist   =  max_rem_dist + 1;

            for (s = dists.begin(); s != dists.end(); s++) {
                if ((*s).second < dist) {
                    min_id = (*s).first;
                    count  = 1;
                    dist   = (*s).second;
                } else if ((*s).second == dist) {
                    count++;
                }
            }

            //
            // Found a merge partner.
            //
            if (min_id >= 0 && count == 1) {
                MergedStack *tag_1 = merged.at(min_id);

                //
                // The max_rem_dist distance allows for the possibility of a frameshift smaller
                // or equal to max_rem_dist at the 3' end of the read. If we detect a possible
                // frameshift, queue the remainder read to re-align the read using the gapped
                // alignment algorithm.
                //
                if (check_frameshift(tag_1, buf, max_rem_dist)) {
                    #pragma omp critical
                    tag_1->rem_queue.push_back(r->id);
                    continue;
                }

                r->utilized = true;
                utilized += r->count();

                #pragma omp critical
                {
                    tag_1->remtags.push_back(r->id);
                    tag_1->count += r->count();
                }
            }
        }

        //
        // Free the k-mers we generated for this remainder read
        //
        for (uint k = 0; k < rem_kmers.size(); k++)
            delete [] rem_kmers[k];

        delete [] buf;
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return utilized;
}

int
search_for_gapped_remainders(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem)
{
    //
    // Search for remainders that can be merged with a gapped alignment.
    //

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    size_t max_rem_len = 0;
    for (auto it = rem.begin(); it != rem.end(); it++) {
        keys.push_back(it->first);
        max_rem_len = it->second->seq->size() > max_rem_len ? it->second->seq->size() : max_rem_len;
    }

    //
    // Calculate the number of k-mers we will generate.
    //
    size_t con_len = 0;
    for (auto miter = merged.begin(); miter != merged.end(); miter++)
        con_len = miter->second->len > con_len ? miter->second->len : con_len;
    size_t kmer_len  = set_kmer_len ? 19 : global_kmer_len;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = (round((double) con_len * min_match_len) - (kmer_len * max_gaps)) - kmer_len + 1;

    //cerr << "  Searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); " << min_hits << " k-mer hits required.\n";

    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);

    #pragma omp parallel
    {
        KmerHashMap::iterator h;
        vector<char *> rem_kmers;
        set<string>    uniq_kmers;
        GappedAln     *aln = new GappedAln(con_len);
        AlignRes       a;
        string         seq, buf;
        char          *rem_buf   = new char[max_rem_len + 1];
        size_t         num_kmers = 0;

        #pragma omp for schedule(dynamic)
        for (uint i = 0; i < keys.size(); i++) {
            auto it = rem.find(keys[i]);
            Rem  *r = it->second;

            if (r->utilized) continue;

            //
            // Generate the k-mers for this remainder sequence
            //
            r->seq->seq(rem_buf);
            if (r->seq->size() < kmer_len) continue;
            num_kmers = r->seq->size() - kmer_len + 1;
            generate_kmers_lazily(rem_buf, kmer_len, num_kmers, rem_kmers);

            uniq_kmers.clear();
            for (uint j = 0; j < num_kmers; j++)
                uniq_kmers.insert(rem_kmers[j]);

            assert(num_kmers > 0);

            map<int, int> hits;
            //
            // Lookup the occurances of each remainder k-mer in the MergedStack k-mer map
            //
            for (set<string>::iterator j = uniq_kmers.begin(); j != uniq_kmers.end(); j++) {
                h = kmer_map.find(j->c_str());

                if (h != kmer_map.end())
                    for (uint k = 0; k <  h->second.size(); k++)
                        hits[h->second[k]]++;
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            vector<Aln> rem_alns;
            Cigar       cigar;

            for (auto hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {

                if (hit_it->second < min_hits) continue;

                MergedStack *tag_1 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_1->masked) continue;

                //
                // Calculate the alignment.
                //
                aln->init(r->seq->size(), tag_1->len);

                if (aln->align(r->seq->seq().c_str(), tag_1->con)) {
                    a = aln->result();
                    rem_alns.push_back(Aln(tag_1->id, a.cigar, a.pct_id, a.gap_cnt));
                }
            }

            if (rem_alns.size() == 0) continue;

            //
            // Rank the alignments for this remainder read.
            //
            sort(rem_alns.begin(), rem_alns.end(), rank_alignments);

            MergedStack *tag_1 = merged[rem_alns[0].id];
            parse_cigar(rem_alns[0].cigar.c_str(), cigar);

            //
            // Check that the alignment still contains fewer than
            // max_utag_dist mismatches.
            //
            if (dist(r->seq->seq().c_str(), tag_1->con, cigar) > max_utag_dist)
                continue;

            //
            // If the alignment has too many gaps, skip it.
            //
            if (rem_alns[0].gap_cnt > (max_gaps + 1))
                continue;

            //
            // If the alignment doesn't span enough of the two sequences, skip it.
            //
            if (rem_alns[0].pct_id < min_match_len)
                continue;

            #pragma omp critical
            tag_1->rem_queue.push_back(r->id);
        }

        //
        // Free the k-mers we generated for this query.
        //
        for (uint j = 0; j < rem_kmers.size(); j++)
            delete [] rem_kmers[j];
        rem_kmers.clear();

        delete [] rem_buf;
        delete aln;
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

size_t
merge_gapped_remainders(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem)
{
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (auto it = merged.begin(); it != merged.end(); it++)
        keys.push_back(it->first);

    size_t utilized = 0;

    #pragma omp parallel reduction(+: utilized)
    {
        Cigar      cigar;
        GappedAln *aln = new GappedAln(merged[keys[0]]->len);
        AlignRes   a;
        string     buf, seq;

        #pragma omp for schedule(dynamic)
        for (uint i = 0; i < keys.size(); i++) {
            MergedStack *tag_1 = merged.at(keys[i]);

            //
            // Don't compute distances for masked tags.
            //
            if (tag_1->masked) continue;

            for (uint i = 0; i < tag_1->rem_queue.size(); i++) {
                Rem *r = rem.at(tag_1->rem_queue[i]);

                // //
                // // We want to align the last max_rem_dist * 2 nucleotides. If sequences are already the wrong length,
                // // adjust so the 3' ends are lined up properly in the alignment matrix.
                // //
                // diff    = tag_1->len - r->seq->size();
                // q_start = r->seq->size() - (max_rem_dist * 2) - 1;
                // q_end   = r->seq->size() - 1;
                // s_start = tag_1->len - (max_rem_dist * 2) - 1;
                // s_end   = tag_1->len - 1;

                // q_start = diff > 0 ? q_start + diff : q_start;
                // s_start = diff < 0 ? s_start + abs(diff) : s_start;

                // cerr << "Consensus size: " << tag_1->len << "; remainder size: " << r->seq->size() << "\n"
                //      << "Con seq: " << tag_1->con << " (" << strlen(tag_1->con) << ")\n"
                //      << "Rem seq: " << r->seq->seq() << "\n"
                //      << "q_start: " << q_start << "; q_end: " << q_end << "; s_start: " << s_start << "; s_end: " << s_end << "\n";

                aln->init(r->seq->size(), tag_1->len);

                if (aln->align(r->seq->seq().c_str(), tag_1->con)) {
                    a = aln->result();
                    parse_cigar(a.cigar.c_str(), cigar);

                    // cerr << "Consensus size: " << tag_1->len << "; remainder size: " << r->seq->size() << "\n"
                    //      << "Con seq: " << tag_1->con << " (" << strlen(tag_1->con) << ")\n"
                    //      << "Rem seq: " << r->seq->seq() << "\n"
                    //      << "Cigar:   " << a.cigar << "\n";

                    //
                    // Check that the alignment still contains fewer than
                    // max_utag_dist mismatches.
                    //
                    if (dist(r->seq->seq().c_str(), tag_1->con, cigar) > max_utag_dist)
                        continue;

                    //
                    // If the alignment has too many gaps, skip it.
                    //
	            if (a.gap_cnt > (max_gaps + 1))
                        continue;

                    //
                    // If the alignment doesn't span enough of the two sequences, skip it.
                    //
                    if (a.pct_id < min_match_len)
                        continue;

                    r->utilized   = true;
                    tag_1->count += r->count();
                    utilized     += r->count();

                    buf = r->seq->seq();
                    seq = apply_cigar_to_seq(buf.c_str(), cigar);
                    r->add_seq(seq.c_str());
                    invert_cigar(cigar);

                    //
                    // Record the gaps (but not soft-masked 3' regions.
                    //
                    uint pos = 0;
                    for (uint j = 0; j < cigar.size(); j++) {
                        if (cigar[j].first == 'I' || cigar[j].first == 'D') {
                            tag_1->gaps.push_back(Gap(pos, pos + cigar[j].second - 1));
                        }
                        pos += cigar[j].second;
                    }
                    edit_gapped_seqs(unique, rem, tag_1, cigar);

                    //
                    // Add this aligned remainder tag to the locus after adjusting the other reads but before adjusting the consensus.
                    //
                    tag_1->remtags.push_back(r->id);
                    update_consensus(tag_1, unique, rem);
                }
            }

            tag_1->rem_queue.clear();
        }
    }

    return utilized;
}

int
call_alleles(MergedStack *mtag, vector<DNANSeq *> &reads, vector<read_type> &read_types)
{
    int     row;
    int     height = reads.size();
    string  allele;
    DNANSeq *d;
    char    base;
    vector<SNP *>::iterator snp;

    for (row = 0; row < height; row++) {
        allele.clear();

        uint snp_cnt = 0;

        //
        // Only call a haplotype from primary reads.
        //
        if (!call_sec_hapl && read_types[row] == secondary) continue;

        for (snp = mtag->snps.begin(); snp != mtag->snps.end(); snp++) {
            if ((*snp)->type != snp_type_het) continue;

            snp_cnt++;

            d    = reads[row];
            base = (*d)[(*snp)->col];

            //
            // Check to make sure the nucleotide at the location of this SNP is
            // of one of the two possible states the multinomial model called.
            //
            if (base == (*snp)->rank_1 || base == (*snp)->rank_2)
                allele += base;
            else
                break;
        }

        if (snp_cnt > 0 && allele.length() == snp_cnt)
            mtag->alleles[allele]++;
    }

    return 0;
}

int
call_consensus(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem, bool invoke_model)
{
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    map<int, MergedStack *>::iterator it;
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++)
        keys.push_back(it->first);

    int i;
    #pragma omp parallel private(i)
    {
        MergedStack *mtag;
        Stack       *utag;
        Rem         *r;

        #pragma omp for schedule(dynamic)
        for (i = 0; i < (int) keys.size(); i++) {
            mtag = merged[keys[i]];

            //
            // Create a two-dimensional array, each row containing one read. For
            // each unique tag that has been merged together, add the sequence for
            // that tag into our array as many times as it originally occurred.
            //
            vector<int>::iterator j;
            vector<DNANSeq *> reads;
            vector<read_type> read_types;
            uint length = unique[mtag->utags.front()]->seq->size();

            for (j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
                utag = unique[*j];

                for (uint k = 0; k < utag->count(); k++) {
                    reads.push_back(utag->seq);
                    read_types.push_back(primary);

                    assert(utag->seq->size() == length);
                }
            }

            // For each remainder tag that has been merged into this Stack, add the sequence.
            for (j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
                r = rem[*j];

                for (uint k = 0; k < r->count(); k++) {
                    reads.push_back(r->seq);
                    read_types.push_back(secondary);

                    assert(r->seq->size() == length);
                }
            }

            //
            // Iterate over each column of the array and call the consensus base.
            //
            uint row, col;
            uint height = reads.size();
            string con;
            map<char, int> nuc;
            map<char, int>::iterator max;
            DNANSeq *d;

            uint cur_gap = mtag->gaps.size() > 0 ? 0 : 1;

            for (col = 0; col < length; col++) {
                nuc['A'] = 0;
                nuc['G'] = 0;
                nuc['C'] = 0;
                nuc['T'] = 0;

                for (row = 0; row < height; row++) {
                    d = reads[row];
                    if (nuc.count((*d)[col]))
                        nuc[(*d)[col]]++;
                }

                //
                // Find the base with a plurality of occurances and call it.
                //
                max = nuc.end();

                for (auto n = nuc.begin(); n != nuc.end(); n++) {
                    if (max == nuc.end() || n->second > max->second)
                        max = n;
                }
                con += max->first;

                //
                // Search this column for the presence of a SNP
                //
                if (invoke_model) {
                    //
                    // Don't invoke the model within gaps.
                    //
                    if (cur_gap < mtag->gaps.size() &&
                        col >= mtag->gaps[cur_gap].start && col <= mtag->gaps[cur_gap].end) {
                        SNP *snp    = new SNP;
                        snp->type   = snp_type_unk;
                        snp->col    = col;
                        snp->rank_1 = '-';
                        snp->rank_2 = '-';
                        mtag->snps.push_back(snp);

                        if (col == mtag->gaps[cur_gap].end)
                            cur_gap++;
                        continue;
                    }

                    switch(model_type) {
                    case snp:
                        call_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case bounded:
                        call_bounded_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case ::fixed:
                        call_multinomial_fixed(mtag, col, nuc);
                        break;
                    default:
                        DOES_NOT_HAPPEN;
                        break;
                    }
                }
            }

            assert(con.length() == length);

            if (invoke_model) {
                assert(mtag->snps.size() == length);

                call_alleles(mtag, reads, read_types);

                if (model_type == ::fixed) {
                    //
                    // Mask nucleotides that are not fixed.
                    //
                    vector<SNP *>::iterator s;
                    for (s = mtag->snps.begin(); s != mtag->snps.end(); s++) {
                        if ((*s)->type == snp_type_unk)
                            con.replace((*s)->col, 1, "N");
                    }
                }
            }

            mtag->add_consensus(con.c_str());
        }
    }

    return 0;
}

int
update_consensus(MergedStack *mtag, map<int, Stack *> &unique, map<int, Rem *> &rem)
{
    Stack *utag;
    Rem   *r;
    //
    // Create a two-dimensional array, each row containing one read. For
    // each unique tag that has been merged together, add the sequence for
    // that tag into our array as many times as it originally occurred.
    //
    vector<DNANSeq *> reads;

    for (auto j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
        utag = unique[*j];

        for (uint k = 0; k < utag->count(); k++)
            reads.push_back(utag->seq);
    }

    // For each remainder tag that has been merged into this Stack, add the sequence.
    for (auto j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
        r = rem[*j];

        for (uint k = 0; k < r->count(); k++)
            reads.push_back(r->seq);
    }

    //
    // Iterate over each column of the array and call the consensus base.
    //
    string   con;
    uint     row, col;
    uint     length = reads[0]->size();
    uint     height = reads.size();
    DNANSeq *d;

    vector<uint> nucs = {0,  // A
                         0,  // C
                         0,  // G
                         0,  // T
                         0}; // N

    con.reserve(length);

    for (col = 0; col < length; col++) {
        nucs.assign(5, 0);

        for (row = 0; row < height; row++) {
            d = reads[row];
            switch ((*d)[col]) {
            case 'A':
                nucs[0]++;
                break;
            case 'C':
                nucs[1]++;
                break;
            case 'G':
                nucs[2]++;
                break;
            case 'T':
                nucs[3]++;
                break;
            case 'N':
                nucs[4]++;
                break;
            }
        }

        //
        // Call the base with a plurality of occurances. Only call 'N' if there is no other information.
        //
        uint max = 0;
        char nuc = 'N';

        for (uint i = 0; i < 4; i++) {
            switch (i) {
            case 0:
                if (nucs[i] > max) {
                    max = nucs[i];
                    nuc = 'A';
                }
                break;
            case 1:
                if (nucs[i] > max) {
                    max = nucs[i];
                    nuc = 'C';
                }
                break;
            case 2:
                if (nucs[i] > max) {
                    max = nucs[i];
                    nuc = 'G';
                }
                break;
            case 3:
                if (nucs[i] > max) {
                    max = nucs[i];
                    nuc = 'T';
                }
                break;
            }
        }
        con += nuc;
    }

    mtag->add_consensus(con.c_str());

    return 0;
}

int
populate_merged_tags(map<int, Stack *> &unique, map<int, MergedStack *> &merged)
{
    map<int, Stack *>::iterator i;
    map<int, MergedStack *>::iterator it_new, it_old;
    Stack       *utag;
    MergedStack *mtag;
    int          k = 0;

    it_old = merged.begin();

    for (i = unique.begin(); i != unique.end(); i++) {
        utag = (*i).second;
        mtag = new MergedStack;

        mtag->id    = k;
        mtag->count = utag->count();
        mtag->utags.push_back(utag->id);
        mtag->add_consensus(utag->seq);

        // Insert the new MergedStack giving a hint as to which position
        // to insert it at.
        it_new = merged.insert(it_old, pair<int, MergedStack *>(k, mtag));
        it_old = it_new;
        k++;
    }

    return 0;
}

void
merge_stacks(map<int, MergedStack *> &merged, size_t& blist_cnt)
{
    map<int, MergedStack *> new_merged;
    map<int, MergedStack *>::const_iterator it;
    map<int, MergedStack *>::iterator it_old, it_new;
    MergedStack *tag_1, *tag_2;
    set<int> merge_map;
    vector<set<int> > merge_lists;

    uint index     = 0;
    int  cohort_id  = 0;
    int  id         = 1;
    uint delev_cnt = 0;
    blist_cnt = 0;

    for (it = merged.begin(); it != merged.end(); it++) {
        tag_1 = it->second;

        //
        // This tag may already have been merged by an earlier operation.
        //
        if (merge_map.count(tag_1->id) > 0)
            continue;

        queue<int>                        merge_list;
        pair<set<int>::iterator,bool>     ret;
        vector<pair<int, int> >::iterator k;

        merge_lists.push_back(set<int>());

        if (tag_1->masked) {
            merge_lists[index].insert(tag_1->id);
            index++;
            continue;
        }

        //
        // Construct a list of MergedStacks that are within a particular distance
        // of this tag.
        //
        merge_lists[index].insert(tag_1->id);
        merge_list.push(tag_1->id);

        while (!merge_list.empty()) {
            tag_2 = merged[merge_list.front()];
            merge_list.pop();

            for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
                ret = merge_lists[index].insert(k->first);

                //
                // If this Tag has not already been added to the merge list (i.e. we were able
                // to insert it in to our unique_merge_list, which is a set), add it for consideration
                // later in the loop.
                //
                if (ret.second == true)
                    merge_list.push((*k).first);
            }
        }

        //
        // Record the nodes that have been merged in this round.
        //
        set<int>::iterator j;
        for (j = merge_lists[index].begin(); j != merge_lists[index].end(); j++)
            merge_map.insert(*j);

        index++;
    }

    #pragma omp parallel private(tag_1, tag_2)
    {
        vector<MergedStack *> merged_tags;

        #pragma omp for reduction(+:delev_cnt) reduction(+:blist_cnt)
        for (uint index = 0; index < merge_lists.size(); index++) {
            //
            // Deal with the simple case of a single locus that does not need to be merged.
            //
            if (merge_lists[index].size() == 1) {
                tag_1 = merged[*(merge_lists[index].begin())];
                tag_2 = merge_tags(merged, merge_lists[index], 0);

                //
                // If this tag is masked, keep the old cohort_id.
                //
                if (tag_1->masked) {
                    tag_2->cohort_id = tag_1->cohort_id;
                } else {
                    tag_2->cohort_id = cohort_id;
                    #pragma omp atomic
                    cohort_id++;
                }
                merged_tags.push_back(tag_2);
                continue;
            }

            //
            // Break large loci down by constructing a minimum
            // spanning tree and severing long distance edges.
            //
            if (deleverage_stacks) {
                vector<MergedStack *> tags;
                bool delev;

                deleverage(merged, merge_lists[index], cohort_id, tags);

                if (tags.size() == 1) {
                    delev = false;
                } else {
                    delev_cnt++;
                    delev = true;
                }

                for (uint t = 0; t < tags.size(); t++) {
                    //tags[t]->id = id;
                    tags[t]->deleveraged = delev;

                    if (tags[t]->utags.size() > max_subgraph) {
                        tags[t]->masked      = true;
                        tags[t]->blacklisted = true;
                        blist_cnt++;
                    }

                    //new_merged.insert(pair<int, MergedStack *>(id, tags[t]));
                    merged_tags.push_back(tags[t]);
                    //id++;
                }

                #pragma omp atomic
                cohort_id++;

            } else {
                //
                // If not deleveraging, merge these tags together into a new MergedStack object.
                //
                tag_2 = merge_tags(merged, merge_lists[index], 0);
                tag_2->cohort_id = cohort_id;

                if (tag_2->utags.size() > max_subgraph) {
                    tag_2->masked      = true;
                    tag_2->blacklisted = true;
                    blist_cnt++;
                }

                //new_merged.insert(pair<int, MergedStack *>(id, tag_2));
                merged_tags.push_back(tag_2);

                #pragma omp atomic
                cohort_id++;
                //id++;
            }
        }

        //
        // Merge the accumulated tags into the new_merged map.
        //
        #pragma omp critical
        {
            it_old = merged.begin();
            for (uint j = 0; j < merged_tags.size(); j++) {
                merged_tags[j]->id = id;
                it_new = new_merged.insert(it_old, pair<int, MergedStack *>(id, merged_tags[j]));
                it_old = it_new;
                id++;
            }
        }
    }

    if (new_merged.empty()) {
        cerr << "Error: Couldn't assemble any loci.\n";
        throw exception();
    }

    //
    // Free the memory from the old map of merged tags.
    //
    for (it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    swap(merged, new_merged);
}

MergedStack *
merge_tags(MergedStack *tag_1, MergedStack *tag_2, int id)
{
    MergedStack *new_tag;

    new_tag     = new MergedStack;
    new_tag->id = id;

    new_tag->deleveraged      = tag_2->deleveraged      || tag_1->deleveraged;
    new_tag->masked           = tag_2->masked           || tag_1->masked;
    new_tag->blacklisted      = tag_2->blacklisted      || tag_1->blacklisted;
    new_tag->gappedlumberjack = tag_2->gappedlumberjack || tag_1->gappedlumberjack;
    new_tag->lumberjackstack  = tag_2->lumberjackstack  || tag_1->lumberjackstack;

    for (uint i = 0; i < tag_1->utags.size(); i++)
        new_tag->utags.push_back(tag_1->utags[i]);

    for (uint i = 0; i < tag_1->remtags.size(); i++)
        new_tag->remtags.push_back(tag_1->remtags[i]);

    for (uint i = 0; i < tag_1->gaps.size(); i++)
        new_tag->gaps.push_back(tag_1->gaps[i]);

    new_tag->count = tag_1->count;

    for (uint i = 0; i < tag_2->utags.size(); i++)
        new_tag->utags.push_back(tag_2->utags[i]);

    for (uint i = 0; i < tag_2->remtags.size(); i++)
        new_tag->remtags.push_back(tag_2->remtags[i]);

    for (uint i = 0; i < tag_2->gaps.size(); i++)
        new_tag->gaps.push_back(tag_2->gaps[i]);

    new_tag->count += tag_2->count;

    return new_tag;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, set<int> &merge_list, int id) {
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (auto i = merge_list.begin(); i != merge_list.end(); i++) {
        tag_2 = merged[(*i)];

        tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
        tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
        tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
        tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

        for (auto j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
            tag_1->utags.push_back(*j);

        for (auto j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
            tag_1->remtags.push_back(*j);

        for (auto j = tag_2->gaps.begin(); j != tag_2->gaps.end(); j++)
            tag_1->gaps.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, int *merge_list, int merge_list_size, int id) {
    int                   i;
    vector<int>::iterator j;
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (i = 0; i < merge_list_size; i++) {
        tag_2 = merged[merge_list[i]];

        tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
        tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
        tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
        tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

        for (j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
            tag_1->utags.push_back(*j);

        for (j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
            tag_1->remtags.push_back(*j);

        for (auto j = tag_2->gaps.begin(); j != tag_2->gaps.end(); j++)
            tag_1->gaps.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

size_t
remove_repetitive_stacks(map<int, MergedStack *> &merged)
{
    //
    // If enabled, check the depth of coverage of each unique tag, and remove
    // from consideration any tags with depths greater than removal_trigger. These tags
    // are likely to be multiple repetitive sites that have been merged together.
    // Because large stacks of unique tags are likely to also generate many one-off
    // sequencing error reads, remove all seqeunces that are a distance of one from
    // the RAD-Tag with high depth of coverage.
    //
    map<int, MergedStack *>::iterator i;
    vector<pair<int, int> >::iterator k;
    map<int, MergedStack *> new_merged;
    MergedStack *tag_1, *tag_2;
    set<int> already_merged;

    const size_t initial_stacks = merged.size();

    //
    // First, iterate through the stacks and populate a list of tags that will be removed
    // (those above the removal_trigger and those 1 nucleotide away). Sort the list of
    // stacks so that we process them from largest depth to shortest so the same stacks
    // are always merged/removed..
    //
    vector<pair<int, int> > ordered_tags;
    for (i = merged.begin(); i != merged.end(); i++) {
        if (i->second->count > removal_trigger)
            ordered_tags.push_back(make_pair(i->second->id, i->second->count));
    }
    sort(ordered_tags.begin(), ordered_tags.end(), compare_pair_intint);

    pair<set<int>::iterator,bool> ret;
    int id = 0;

    //
    // Merge all stacks that are over the removal trigger with their nearest neighbors and
    // mask them so they are not further processed by the program.
    //
    for (uint j = 0; j < ordered_tags.size(); j++) {
        tag_1 = merged[ordered_tags[j].first];

        //
        // Don't process a tag that has already been merged.
        //
        if (already_merged.count(tag_1->id) > 0)
            continue;

        //
        // Construct a list of MergedStacks that are either:
        //   within a distance of 1 nucleotide of this tag, or
        //   are they themselves above the lumberjack stacks limit.
        //
        queue<int> merge_queue;
        set<int>   merge_list;
        merge_queue.push(tag_1->id);
        merge_list.insert(tag_1->id);
        already_merged.insert(tag_1->id);

        while (!merge_queue.empty()) {
            tag_2 = merged[merge_queue.front()];
            merge_queue.pop();

            if (tag_2->count < removal_trigger)
                continue;

            for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
                ret = already_merged.insert(k->first);

                if (ret.second == true) {
                    merge_queue.push(k->first);
                    merge_list.insert(k->first);
                }
            }
        }

        //
        // Merge these tags together into a new MergedStack object.
        //
        tag_2 = merge_tags(merged, merge_list, id);
        tag_2->add_consensus(tag_1->con);

        tag_2->lumberjackstack = true;
        tag_2->masked          = true;
        tag_2->blacklisted     = true;

        new_merged.insert(make_pair(id, tag_2));
        id++;
    }

    //
    // Move the non-lumberjack stacks, unmodified, into the new merged map.
    //
    size_t remaining_stacks = 0;
    for (i = merged.begin(); i != merged.end(); i++) {
        tag_1 = i->second;

        if (already_merged.count(tag_1->id) > 0)
            continue;
        ++remaining_stacks;

        set<int> merge_list;
        merge_list.insert(tag_1->id);

        tag_2 = merge_tags(merged, merge_list, id);
        tag_2->add_consensus(tag_1->con);

        new_merged.insert(make_pair(id, tag_2));
        id++;
    }

    //
    // Free the memory from the old map of merged tags.
    //
    map<int, MergedStack *>::iterator it;
    for (it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    merged.swap(new_merged);
    return initial_stacks - remaining_stacks;
}

int deleverage(map<int, MergedStack *> &merged,
               set<int> &merge_list,
               int cohort_id,
               vector<MergedStack *> &deleveraged_tags) {
    set<int>::iterator it;
    vector<pair<int, int> >::iterator j;
    MergedStack *tag_1, *tag_2;
    uint k, l;

    //
    // Create a minimum spanning tree in order to determine the minimum distance
    // between each node in the list.
    //
    mst::MinSpanTree *mst = new mst::MinSpanTree();
    vector<int>  keys;

    for (it = merge_list.begin(); it != merge_list.end(); it++) {
        keys.push_back(*it);
        mst->add_node(*it);
        tag_1 = merged[*it];

        // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    //
    // Measure the distance between each pair of nodes and add edges to our
    // minimum spanning tree.
    //
    mst::Node *n_1, *n_2;
    for (k = 0; k < keys.size(); k++) {
        tag_1 = merged[keys[k]];
        n_1   = mst->node(keys[k]);

        for (l = k+1; l < keys.size(); l++) {
            tag_2 = merged[keys[l]];
            n_2   = mst->node(keys[l]);

            int d = dist(tag_1, tag_2);

            n_1->add_edge(mst->node(keys[l]), d);
            n_2->add_edge(mst->node(keys[k]), d);
        }
    }

    //
    // Build the minimum spanning tree.
    //
    mst->build_tree();

    //
    // Visualize the MST
    //
    if (dump_graph) {
        stringstream gout_file;
        gout_file << prefix_path << "_" << keys[0] << ".dot";
        string vis = mst->vis(true);
        ofstream gvis(gout_file.str().c_str());
        gvis << vis;
        gvis.close();
    }

    set<int> visited;
    set<int, int_increasing> dists;
    queue<mst::Node *> q;

    mst::Node *n = mst->head();
    q.push(n);

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);
                // cerr << n->id << " -> " << n->min_adj_list[i]->id << ": ";

                //
                // Find the edge distance.
                //
                for (uint j = 0; j < n->edges.size(); j++)
                    if (n->edges[j]->child == n->min_adj_list[i]) {
                        // cerr << n->edges[j]->dist << "\n";
                        dists.insert(n->edges[j]->dist);
                    }
            }
        }
    }

    //
    // This set is sorted by definition. Check if there is more than a single
    // distance separating stacks.
    //
    if (dists.size() == 1) {
        tag_1 = merge_tags(merged, merge_list, 0);
        deleveraged_tags.push_back(tag_1);
        delete mst;
        return 0;
    }

    uint min_dist = *(dists.begin());

    //
    // If there is more than a single distance, split the minimum spanning tree
    // into separate loci, by cutting the tree at the larger distances.
    //
    set<int> uniq_merge_list;
    visited.clear();
    n = mst->head();
    q.push(n);
    int id = 0;

    uniq_merge_list.insert(n->id);
    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);

                for (uint j = 0; j < n->edges.size(); j++) {
                    if (n->edges[j]->child == n->min_adj_list[i])
                        if (n->edges[j]->dist > min_dist) {

                            // cerr << "Merging the following stacks into a locus:\n";
                            for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
                                tag_1 = merged[*it];
                                // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
                            }

                            tag_1 = merge_tags(merged, uniq_merge_list, id);
                            tag_1->cohort_id = cohort_id;
                            deleveraged_tags.push_back(tag_1);
                            uniq_merge_list.clear();
                            id++;
                        }
                }

                uniq_merge_list.insert(n->min_adj_list[i]->id);
            }
        }
    }

    // cerr << "Merging the following stacks into a locus:\n";
    for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
        tag_1 = merged[*it];
        // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    tag_1 = merge_tags(merged, uniq_merge_list, id);
    tag_1->cohort_id = cohort_id;
    deleveraged_tags.push_back(tag_1);
    uniq_merge_list.clear();

    delete mst;

    return 0;
}

int
calc_kmer_distance(map<int, MergedStack *> &merged, int utag_dist)
{
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    MergedStack   *tag_1, *tag_2;
    map<int, MergedStack *>::iterator it;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++)
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    size_t con_len   = 0;
    for (auto miter = merged.begin(); miter != merged.end(); miter++)
        con_len = miter->second->len > con_len ? miter->second->len : con_len;
    size_t kmer_len  = set_kmer_len ? determine_kmer_length(con_len, utag_dist) : global_kmer_len;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, utag_dist, con_len, set_kmer_len ? true : false);

    //cerr << "  Distance allowed between stacks: " << utag_dist
    //     << "; searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); "
    //     << min_hits << " k-mer hits required.\n";

    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);

    #pragma omp parallel private(tag_1, tag_2)
    {
        KmerHashMap::iterator h;
        vector<char *>        query_kmers;

        size_t num_kmers = con_len - kmer_len + 1;
        initialize_kmers(kmer_len, num_kmers, query_kmers);

        #pragma omp for schedule(dynamic)
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            if (tag_1->len < kmer_len) continue;
            num_kmers = tag_1->len - kmer_len + 1;
            generate_kmers_lazily(tag_1->con, kmer_len, num_kmers, query_kmers);

            map<int, int> hits;
            int d;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (uint j = 0; j < num_kmers; j++) {
                h = kmer_map.find(query_kmers[j]);

                if (h != kmer_map.end())
                    for (uint k = 0; k <  h->second.size(); k++)
                        hits[h->second[k]]++;
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {

                if (hit_it->second < min_hits) continue;

                tag_2 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                d = dist(tag_1, tag_2);

                //
                // Check if any of the mismatches occur at the 3' end of the read. If they
                // do, they may indicate a frameshift is present at the 3' end of the read,
                // which will cause problems when we try to merge loci across samples.
                // If found, do not merge these tags, leave them for the gapped alignmnet
                // algorithm.
                //
                if (d <= utag_dist && check_frameshift(tag_1, tag_2, (size_t) utag_dist))
                    continue;

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance. Thesesequences will then be
                //  merged in the following step of the algorithm.
                //
                if (d <= utag_dist)
                    tag_1->add_dist(tag_2->id, d);
            }
            // Sort the vector of distances.
            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
        }

        //
        // Free the k-mers we generated for this query
        //
        for (uint j = 0; j < query_kmers.size(); j++)
            delete [] query_kmers[j];
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

//int calc_distance(map<int, MergedStack *> &merged, int utag_dist) {
//    //
//    // Calculate the distance (number of mismatches) between each pair
//    // of Radtags. We expect all radtags to be the same length;
//    //
//    map<int, MergedStack *>::iterator it;
//    MergedStack *tag_1, *tag_2;
//    int i, j;
//
//    // OpenMP can't parallelize random access iterators, so we convert
//    // our map to a vector of integer keys.
//    vector<int> keys;
//    for (it = merged.begin(); it != merged.end(); it++)
//        keys.push_back(it->first);
//
//    #pragma omp parallel private(i, j, tag_1, tag_2)
//    {
//        #pragma omp for schedule(dynamic)
//        for (i = 0; i < (int) keys.size(); i++) {
//
//            tag_1 = merged[keys[i]];
//
//            // Don't compute distances for masked tags
//            if (tag_1->masked) continue;
//
//            int d;
//
//            for (j = 0; j < (int) keys.size(); j++) {
//                tag_2 = merged[keys[j]];
//
//                // Don't compute distances for masked tags
//                if (tag_2->masked) continue;
//
//                // Don't compare tag_1 against itself.
//                if (tag_1 == tag_2) continue;
//
//                d = dist(tag_1, tag_2);
//                //cerr << "  Distance: " << d << "\n";
//
//                //
//                // Store the distance between these two sequences if it is
//                // below the maximum distance (which governs those
//                // sequences to be merged in the following step of the
//                // algorithm.)
//                //
//                if (d == utag_dist) {
//                    tag_1->add_dist(tag_2->id, d);
//                    //cerr << "  HIT.\n";
//                }
//            }
//
//            // Sort the vector of distances.
//            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
//        }
//    }
//
//    return 0;
//}

void reduce_radtags(
        DNASeqHashMap &radtags, map<int, Stack *> &unique,
        map<int, Rem *> &rem,
        size_t& n_u_reads,
        size_t& n_r_reads
        ) {
    n_u_reads = 0;
    n_r_reads = 0;
    DNASeqHashMap::iterator it;

    Rem   *r;
    Stack *u;
    int   global_id = 1;

    for (it = radtags.begin(); it != radtags.end(); it++) {
        if (it->second.count() < min_merge_cov) {
            //
            // Don't record this unique RAD-Tag if its coverage is below
            // the specified cutoff. However, add the reads to the remainder
            // vector for later processing.
            //
            r     = new Rem;
            r->id = global_id;
            r->add_seq(&it->first);

            for (uint i = 0; i < it->second.ids.size(); i++)
                r->add_id(it->second.ids[i]);
            n_r_reads += r->map.size();

            rem[r->id] = r;

        } else {
            //
            // Populate a Stack object for this unique radtag. Create a
            // map of the IDs for the sequences that have been
            // collapsed into this radtag.
            //
            u     = new Stack;
            u->id = global_id;
            u->add_seq(&it->first);

            // Copy the original Fastq IDs from which this unique radtag was built.
            for (uint i = 0; i < it->second.ids.size(); i++)
                u->add_id(it->second.ids[i]);
            n_u_reads += u->map.size();

            unique[u->id] = u;
        }
        if (global_id == INT_MAX)
            throw std::overflow_error("overflow in reduce_radtags()");
        ++global_id;
    }

    if (unique.size() == 0) {
        cerr << "Error: Unable to form any primary stacks.\n";
        exit(1);
    }
}

void
calc_coverage_distribution(map<int, MergedStack *> &merged,
                           double &mean, double &stdev, double &max, double &tot_depth)
{
    max   = 0.0;
    mean = 0.0;
    stdev = 0.0;
    tot_depth = 0.0;

    size_t not_blacklisted = 0;
    for (const pair<int, MergedStack*>& mtag : merged) {
        if (mtag.second->blacklisted)
            continue;
        ++not_blacklisted;

        double depth = mtag.second->count;
        tot_depth += depth;
        if (depth > max)
            max = depth;
    }
    mean = tot_depth / not_blacklisted;

    for (const pair<int, MergedStack*>& mtag : merged)
        if (!mtag.second->blacklisted)
            stdev += pow(mtag.second->count - mean, 2);
    stdev /= not_blacklisted;
    stdev = sqrt(stdev);
}

int
write_results(map<int, MergedStack *> &m, map<int, Stack *> &u, map<int, Rem *> &r)
{
    map<int, MergedStack *>::iterator i;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack  *tag_1;
    Stack        *tag_2;
    Rem          *rem;
    char strbuf32[32];

    bool gzip = (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) ? true : false;

    //
    // Read in the set of sequencing IDs so they can be included in the output.
    //
    vector<char *> seq_ids;
    load_seq_ids(seq_ids);

    //
    // Create the output file names.
    //
    string tag_file = prefix_path + ".tags.tsv";
    string snp_file = prefix_path + ".snps.tsv";
    string all_file = prefix_path + ".alleles.tsv";
    if (gzip) {
        tag_file += ".gz";
        snp_file += ".gz";
        all_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    VersatileWriter tags{tag_file}, snps{snp_file}, alle{all_file};

    //
    // Record the version of Stacks used and the date generated as a comment in the catalog.
    //
    // Obtain the current date.
    //
    stringstream log;
    time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log << "# ustacks version " << VERSION << "; generated on " << date << "\n";
    tags << log.str();
    snps << log.str();
    alle << log.str();

    for (i = m.begin(); i != m.end(); i++) {
        float total = 0;
        tag_1 = i->second;

        //
        // Calculate the log likelihood of this merged stack.
        //
        tag_1->gen_matrix(u, r);
        tag_1->calc_likelihood();

        // First write the consensus sequence
        tags << sample_id          << "\t"
             << tag_1->id          << "\t"
             << "consensus\t"      << "\t"
             << "\t"
             << tag_1->con         << "\t"
             << tag_1->deleveraged << "\t"
             << tag_1->blacklisted << "\t"
             << tag_1->lumberjackstack << "\n";

        //
        // Write a sequence recording the output of the SNP model for each nucleotide.
        //
        tags << sample_id << "\t"
             << tag_1->id << "\t"
             << "model\t" << "\t"
             << "\t";
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            switch((*s)->type) {
            case snp_type_het:
                tags << "E";
                break;
            case snp_type_hom:
                tags << "O";
                break;
            default:
                tags << "U";
                break;
            }
        }
        tags << "\t\t\t\n";

        //
        // Now write out the components of each unique tag merged into this locus.
        //
        int id = 0;
        for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++) {
            tag_2  = u[*k];
            total += tag_2->count();

            for (uint j = 0; j < tag_2->map.size(); j++) {
                tags << sample_id << "\t"
                     << tag_1->id << "\t"
                     << "primary\t"
                     << id << "\t"
                     << seq_ids[tag_2->map[j]] << "\t"
                     << tag_2->seq->seq()
                     << "\t\t\t\n";
            }
            id++;
        }

        //
        // Write out the remainder tags merged into this unique tag.
        //
        for (k = tag_1->remtags.begin(); k != tag_1->remtags.end(); k++) {
            rem    = r[*k];
            total += rem->map.size();
            for (uint j = 0; j < rem->map.size(); j++)
                tags << sample_id << "\t"
                     << tag_1->id << "\t"
                     << "secondary\t"
                     << "\t"
                     << seq_ids[rem->map[j]] << "\t"
                     << rem->seq->seq()
                     << "\t\t\t\n";
        }

        //
        // Write out the model calls for each nucleotide in this locus.
        //
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            snps << sample_id << "\t"
                 << tag_1->id << "\t"
                 << (size_t) (*s)->col << "\t";
            switch((*s)->type) {
            case snp_type_het:
                snps << "E\t";
                break;
            case snp_type_hom:
                snps << "O\t";
                break;
            default:
                snps << "U\t";
                break;
            }
            sprintf(strbuf32, "%.2f", (*s)->lratio);
            snps << strbuf32 << "\t"
                 << (*s)->rank_1 << "\t"
                 << (*s)->rank_2 << "\t\t\n";
        }

        //
        // Write the expressed alleles seen for the recorded SNPs and
        // the percentage of tags a particular allele occupies.
        //
        for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            alle << sample_id   << "\t"
                 << tag_1->id   << "\t"
                 << (*t).first  << "\t";
            sprintf(strbuf32, "%.2f", ((*t).second/total) * 100);
            alle << strbuf32 << "\t"
                 << (*t).second << "\n";
        }
    }

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log.str("");
    log << "# ustacks completed on " << date << "\n";
    tags << log.str();
    snps << log.str();
    alle << log.str();

    tags.close();
    snps.close();
    alle.close();

    //
    // Free sequence IDs.
    //
    for (uint i = 0; i < seq_ids.size(); i++)
        delete [] seq_ids[i];

    //
    // If specified, output reads not utilized in any stacks.
    //
    if (retain_rem_reads) {
        string unused_file = prefix_path + ".unused.fa";
        VersatileWriter unused{unused_file};

        map<int, Rem *>::iterator r_it;
        for (r_it = r.begin(); r_it != r.end(); r_it++)
            if (r_it->second->utilized == false)
                unused << ">" << (long) r_it->second->id << "\n" << r_it->second->seq->seq() << "\n";
        unused.close();
    }

    return 0;
}

int dump_stack_graph(string data_file,
                     map<int, Stack *> &unique,
                     map<int, MergedStack *> &merged,
                     vector<int> &keys,
                     map<int, map<int, double> > &dist_map,
                     map<int, set<int> > &cluster_map) {
    uint s, t;
    double d, scale, scaled_d;
    char label[32];
    vector<string> colors;
    ofstream data(data_file.c_str());

    size_t pos_1 = data_file.find_last_of("/");
    size_t pos_2 = data_file.find_last_of(".");
    string title = data_file.substr(pos_1 + 1, pos_2 - pos_1 - 1);

    //
    // Output a list of IDs so we can locate these stacks in the final results.
    //
    for (s = 0; s < keys.size(); s++)
        data << "/* " << keys[s] << ": " << unique[merged[keys[s]]->utags[0]]->map[0] << "; depth: " << merged[keys[s]]->count << " */\n";

    //
    // Output a specification to visualize the stack graph using graphviz:
    //   http://www.graphviz.org/
    //
    data << "graph " << title.c_str() << " {\n"
         << "rankdir=LR\n"
         << "size=\"20!\"\n"
         << "overlap=false\n"
         << "node [shape=circle style=filled fillcolor=\"#3875d7\" fontname=\"Arial\"];\n"
         << "edge [fontsize=8.0 fontname=\"Arial\" color=\"#aaaaaa\"];\n";

    colors.push_back("red");
    colors.push_back("blue");
    colors.push_back("green");
    colors.push_back("brown");
    colors.push_back("purple");

    map<int, set<int> >::iterator c;
    set<int>::iterator it;
    int color_index = 0;
    string color;

    // Write out the clusters created by R, prior to writing all the nodes and connections.
    s = 0;
    for (c = cluster_map.begin(); c != cluster_map.end(); c++) {
        data << "subgraph " << s << " {\n"
             << "    edge [penwidth=5 fontsize=12.0 fontcolor=\"black\" color=\"black\"]\n";

        if ((*c).second.size() == 1) {
            color = "white";
            data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"black\"]\n";
        } else {
            color = colors[color_index % colors.size()];
            data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"white\"]\n";
            color_index++;
        }

        for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
            data << "    " << *it << "\n";
        }

        if ((*c).second.size() > 1) {
            uint j = 0;
            for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
                data << *it;
                if (j < (*c).second.size() - 1)
                    data << " -- ";
                j++;
            }
        }

        data << "}\n";
        s++;
    }

    //
    // Scale the graph to display on a 10 inch canvas. Find the largest edge weight
    // and scale the edge lengths to fit the canvas.
    //
    scale = 0.0;
    for (s = 0; s < keys.size(); s++)
        for (t = s+1; t < keys.size(); t++)
            scale = dist_map[keys[s]][keys[t]] > scale ? dist_map[keys[s]][keys[t]] : scale;
    scale = scale / 20;

    for (s = 0; s < keys.size(); s++) {
        for (t = s+1; t < keys.size(); t++) {
            d = dist_map[keys[s]][keys[t]];
            scaled_d = d / scale;
            scaled_d = scaled_d < 0.75 ? 0.75 : scaled_d;
            sprintf(label, "%.1f", d);
            data << keys[s] << " -- " << keys[t] << " [len=" << scaled_d << ", label=" << label << "];\n";
        }
    }

    data << "}\n";

    data.close();

    return 0;
}

int dump_unique_tags(map<int, Stack *> &u) {
    map<int, Stack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator mit;

    for (it = u.begin(); it != u.end(); it++) {
        cerr << "UniqueTag UID: " << (*it).second->id << "\n"
             << "  Seq:       "   << it->second->seq->seq() << "\n"
             << "  IDs:       ";

        for (uint j = 0; j < it->second->map.size(); j++)
            cerr << it->second->map[j] << " ";

        cerr << "\n\n";
    }

    return 0;
}

int dump_merged_tags(map<int, MergedStack *> &m) {
    map<int, MergedStack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator fit;

    for (it = m.begin(); it != m.end(); it++) {

        cerr << "MergedStack ID: " << it->second->id << "\n"
             << "  Consensus:  ";
        if (it->second->con != NULL)
            cerr << it->second->con << "\n";
        else
            cerr << "\n";
        cerr << "  IDs:        ";

        for (fit = it->second->utags.begin(); fit != it->second->utags.end(); fit++)
            cerr << (*fit) << " ";

        cerr << "\n"
             << "  Distances: ";

        for (pit = it->second->dist.begin(); pit != it->second->dist.end(); pit++)
            cerr << (*pit).first << ": " << (*pit).second << ", ";

        cerr << "\n\n";
    }

    return 0;
}

void load_radtags(string in_file, DNASeqHashMap &radtags, size_t& n_reads) {
    n_reads = 0;

    Input *fh = NULL;

    if (in_file_type == FileT::fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == FileT::fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(in_file.c_str());

    long  int corrected = 0;
    size_t i            = 0;
    short int seql      = 0;
    short int prev_seql = 0;
    bool  len_mismatch  = false;

    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
        if (i % 1000000 == 0 && i>0)
            cerr << i/1000000 << "M..." << flush;

        prev_seql = seql;
        seql      = 0;

        for (char *p = c.seq; *p != '\0'; p++, seql++)
            switch (*p) {
            case 'N':
            case 'n':
            case '.':
                *p = 'A';
                corrected++;
            }

        if (seql != prev_seql && prev_seql > 0)
            len_mismatch = true;

        DNASeqHashMap::iterator element = radtags.insert({DNANSeq(seql, c.seq), HVal()}).first;
        element->second.add_id(i);
        i++;
    }
    cerr << "\n";

    if (i == 0) {
        cerr << "Error: Unable to load data from '" << in_file.c_str() << "'.\n";
        exit(1);
    }
    if (len_mismatch && !force_diff_len) {
        cerr << "Error: different sequence lengths detected, this will interfere with Stacks "
             << "algorithms, trim reads to uniform length (override this check with --force-diff-len).\n";
        exit(1);
    } else if (force_diff_len) {
        cerr << "Warning: different sequence lengths detected, this could interfere with Stacks algorithms.\n";
    }
    if (corrected > 0)
        cerr << "Warning: Input reads contained " << corrected << " uncalled nucleotides.\n";

    delete fh;
    n_reads = i;
}

int
load_seq_ids(vector<char *> &seq_ids)
{
    Input *fh = NULL;

    if (in_file_type == FileT::fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == FileT::fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(in_file.c_str());

    cerr << "Refetching read IDs...";

    char *id;
    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
        id = new char[strlen(c.id) + 1];
        strcpy(id, c.id);
        seq_ids.push_back(id);
    }

    delete fh;
    return 0;
}

long double factorial(int i) {
    long double f = 1;

    if (i == 0) return 1;

    do {
        f = f * i;
        i--;
    } while (i > 0);

    return f;
}

int parse_command_line(int argc, char* argv[]) {
    string out_path;
    string sample_name;

    int c;
    while (1) {
        static struct option long_options[] = {
            {"help",             no_argument,       NULL, 'h'},
            {"version",          no_argument,       NULL, 'v'},
            {"infile-type",      required_argument, NULL, 't'}, {"infile_type",      required_argument, NULL, 't'},
            {"file",             required_argument, NULL, 'f'},
            {"outpath",          required_argument, NULL, 'o'},
            {"name",             required_argument, NULL, 1002},
            {"id",               required_argument, NULL, 'i'},
            {"min-cov",          required_argument, NULL, 'm'}, {"min_cov",          required_argument, NULL, 'm'},
            {"max-dist",         required_argument, NULL, 'M'}, {"max_dist",         required_argument, NULL, 'M'},
            {"max-sec-dist",     required_argument, NULL, 'N'}, {"max_sec_dist",     required_argument, NULL, 'N'},
            {"max-locus-stacks", required_argument, NULL, 'K'}, {"max_locus_stacks", required_argument, NULL, 'K'},
            {"k-len",            required_argument, NULL, 'k'}, {"k_len",            required_argument, NULL, 'k'},
            {"num-threads",      required_argument, NULL, 'p'}, {"num_threads",      required_argument, NULL, 'p'},
            {"deleverage",       no_argument,       NULL, 'd'},
            {"keep-high-cov",    no_argument,       NULL, 1000}, {"keep_high_cov",    no_argument,       NULL, 1000},
            {"high-cov-thres",   required_argument, NULL, 1001}, {"high_cov_thres",   required_argument, NULL, 1001},
            {"retain-rem",       no_argument,       NULL, 'R'}, {"retain_rem",       no_argument,       NULL, 'R'},
            {"graph",            no_argument,       NULL, 'g'},
            {"sec-hapl",         no_argument,       NULL, 'H'}, {"sec_hapl",         no_argument,       NULL, 'H'},
            {"disable-gapped",   no_argument,       NULL, 'G'},
            {"max-gaps",         required_argument, NULL, 'X'}, {"max_gaps",         required_argument, NULL, 'X'},
            {"min-aln-len",      required_argument, NULL, 'x'}, {"min_aln_len",      required_argument, NULL, 'x'},
            {"model-type",       required_argument, NULL, 'T'}, {"model_type",       required_argument, NULL, 'T'},
            {"bc-err-freq",      required_argument, NULL, 'e'}, {"bc_err_freq",      required_argument, NULL, 'e'},
            {"bound-low",        required_argument, NULL, 'L'}, {"bound_low",        required_argument, NULL, 'L'},
            {"bound-high",       required_argument, NULL, 'U'}, {"bound_high",       required_argument, NULL, 'U'},
            {"alpha",            required_argument, NULL, 'A'},
            {"r-deprecated",     no_argument,       NULL, 'r'},
            {"force-diff-len",   no_argument,       NULL, 1003}, {"force_diff_len",   no_argument,       NULL, 1003},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "GhHvdrgRA:L:U:f:o:s:i:m:e:p:t:M:N:K:k:T:X:x:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 't':
            if (strcmp(optarg, "tsv") == 0)
                in_file_type = FileT::tsv;
            else if (strcmp(optarg, "fasta") == 0)
                in_file_type = FileT::fasta;
            else if (strcmp(optarg, "fastq") == 0)
                in_file_type = FileT::fastq;
            else if (strcasecmp(optarg, "gzfasta") == 0)
                in_file_type = FileT::gzfasta;
            else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = FileT::gzfastq;
            else
                in_file_type = FileT::unknown;
            break;
        case 'f':
            in_file = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 1002:
            sample_name = optarg;
            break;
        case 'i':
            sample_id = is_integer(optarg);
            break;
        case 'm':
            min_merge_cov = is_integer(optarg);
            break;
        case 'M':
            max_utag_dist = is_integer(optarg);
            break;
        case 'N':
            max_rem_dist = is_integer(optarg);
            break;
        case 'd':
            deleverage_stacks++;
            break;
        case 1000: // keep_high_cov
            remove_rep_stacks = false;
            break;
        case 1001: // high_cov_thres
            removal_threshold = is_double(optarg);
            if (removal_threshold <= 0.0) {
                cerr << "Error: Bad threshold: '" << optarg << "'\n";
                help();
            }
            break;
        case 'K':
            max_subgraph = is_integer(optarg);
            break;
        case 'k':
            set_kmer_len    = false;
            global_kmer_len = is_integer(optarg);
            break;
        case 'R':
            retain_rem_reads = true;
            break;
        case 'g':
            dump_graph++;
            break;
        case 'G':
            gapped_alignments = false;
            break;
        case 'X':
            max_gaps = is_double(optarg);
            break;
        case 'x':
            min_match_len = is_double(optarg);
            break;
        case 'T':
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = ::fixed;
            } else if (strcmp(optarg, "bounded") == 0) {
                model_type = bounded;
            } else {
                cerr << "Unknown model type specified '" << optarg << "'\n";
                help();
            }
            break;;
        case 'e':
            barcode_err_freq = is_double(optarg);
            break;
        case 'L':
            bound_low  = is_double(optarg);
            break;
        case 'U':
            bound_high = is_double(optarg);
            break;
        case 'A':
            alpha = is_double(optarg);
            break;
        case 'H':
            call_sec_hapl = false;
            break;
        case 'p':
            num_threads = is_integer(optarg);
            break;
        case 'v':
            version();
            break;
        case 1003:
            force_diff_len = true;
            break;
        case 'r': // deprecated Dec 2016, v1.45
            cerr << "Warning: Ignoring deprecated option -r (this has become the default).\n";
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;

        default:
            cerr << "Unknown command line option '" << (char) c << "'\n";
            help();
            exit(1);
        }
    }

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    if (set_kmer_len == false && (global_kmer_len < 7 || global_kmer_len > 31)) {
        cerr << "Kmer length must be between 7 and 31bp.\n";
        help();
    }

    if (alpha != 0.1 && alpha != 0.05 && alpha != 0.01 && alpha != 0.001) {
        cerr << "SNP model alpha significance level must be either 0.1, 0.05, 0.01, or 0.001.\n";
        help();
    }

    if (bound_low != 0 && (bound_low < 0 || bound_low >= 1.0)) {
        cerr << "SNP model lower bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_high != 1 && (bound_high <= 0 || bound_high > 1.0)) {
        cerr << "SNP model upper bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_low > 0 || bound_high < 1.0) {
        model_type = bounded;
    }

    if (in_file.empty()) {
        cerr << "You must specify an input file.\n";
        help();
    }

    if (in_file_type == FileT::unknown) {
        in_file_type = guess_file_type(in_file);
        if (in_file_type == FileT::unknown) {
            cerr << "Unable to recongnize the extention of file '" << in_file << "'.\n";
            help();
        }
    }

    // Set `prefix_path`.
    if (out_path.length() == 0)
        out_path = ".";
    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";
    if (sample_name.empty()) {
        size_t pos_1 = in_file.find_last_of("/");
        size_t pos_2 = in_file.find_last_of(".");
        if (in_file.substr(pos_2) == ".gz")
            pos_2 = in_file.substr(0, pos_2).find_last_of(".");
        if (pos_2 == string::npos || pos_2 <= pos_1 + 1) {
            cerr << "Failed to guess the sample's name.\n";
            help();
        }
        sample_name = in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1));
    }
    prefix_path = out_path + sample_name;

    if (sample_id < 0) {
        cerr << "A unique sample ID must be provided (a positive integer).\n";
        help();
    }

    if (model_type == ::fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        help();
    }

    return 0;
}

void version() {
    cerr << "ustacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "ustacks " << VERSION << "\n"
         << "ustacks -f file_path -i id -o path [-M max_dist] [-m min_cov] [-p num_threads]" << "\n"
         << "  f: input file path.\n"
         << "  i: a unique integer ID for this sample.\n"
         << "  o: output path to write results.\n"
         << "  M: Maximum distance (in nucleotides) allowed between stacks (default 2).\n"
         << "  m: Minimum depth of coverage required to create a stack (default 3).\n"
         << "  N: Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).\n"
         << "  p: enable parallel execution with num_threads threads.\n"
         << "  t: input file type. Supported types: fasta, fastq, gzfasta, or gzfastq (default: guess).\n"
         << "  --name: a name for the sample (default: input file name minus the suffix).\n"
         << "  R: retain unused reads.\n"
         << "  H: disable calling haplotypes from secondary reads.\n"
         << "\n"
         << "  Stack assembly options:\n"
         << "    d,--deleverage: enable the Deleveraging algorithm, used for resolving over merged tags.\n"
         << "    --keep-high-cov: disable the algorithm that removes highly-repetitive stacks and nearby errors.\n"
         << "    --high-cov-thres: highly-repetitive stacks threshold, in standard deviation units (default: 3.0).\n"
         << "    --max-locus-stacks <num>: maximum number of stacks at a single de novo locus (default 3).\n"
         << "     --k-len <len>: specify k-mer size for matching between alleles and loci (automatically calculated by default).\n\n"
         << "  Gapped assembly options:\n"
         << "    --max-gaps: number of gaps allowed between stacks before merging (default: 2).\n"
         << "    --min-aln-len: minimum length of aligned sequence in a gapped alignment (default: 0.80).\n\n"
         << "    --disable-gapped: do not preform gapped alignments between stacks (default: gapped alignements enabled).\n"
         << "  Model options:\n"
         << "    --model-type: either 'snp' (default), 'bounded', or 'fixed'\n"
         << "    For the SNP or Bounded SNP model:\n"
         << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
         << "    For the Bounded SNP model:\n"
         << "      --bound-low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
         << "      --bound-high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
         << "    For the Fixed model:\n"
         << "      --bc-err-freq <num>: specify the barcode error frequency, between 0 and 1.0.\n"
         << "\n"
         << "  h: display this help messsage.\n";

    exit(1);
}
