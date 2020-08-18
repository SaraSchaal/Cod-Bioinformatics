// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2015, Julian Catchen <jcatchen@illinois.edu>
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
// phasedstacks -- analyse phased data, descended from a Stacks analysis.
//

#include "phasedstacks.h"

// Global variables to hold command-line options.
FileT  in_file_type = FileT::unknown;
int    num_threads  = 1;
int    batch_id     = 0;
string cat_path;
string in_path;
string out_path;
string out_file;
string pmap_path;
bool   haplotypes       = false;
bool   write_zeros      = true;
double p_value_cutoff   = 0.05;
double chi_sq_limit     = 3.84;
double minor_freq_lim   = 0.1;
double min_inform_pairs = 0.90;
uint   max_pair_dist    = 1000000;
uint   bucket_dist      = 5000;
double dprime_threshold = false;
double dprime_threshold_level = 0.0;

set<int> whitelist, blacklist;

map<string, int> pop_map;
map<int, int>    pop_cnts;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    if (p_value_cutoff == 0.1) {
        chi_sq_limit =  2.71;
    } else if (p_value_cutoff == 0.05) {
        chi_sq_limit =  3.84;
    } else if (p_value_cutoff == 0.01) {
        chi_sq_limit =  6.64;
    } else if (p_value_cutoff == 0.001) {
        chi_sq_limit = 10.83;
    }

    cerr << "Minor allele frequency cutoff: " << minor_freq_lim << "\n"
         << "Looking for ";
    switch(in_file_type) {
    case FileT::beagle:
        cerr << "Beagle";
        break;
    case FileT::phase:
        cerr << "PHASE";
        break;
    case FileT::fastphase:
    default:
        cerr << "fastPhase";
        break;
    }
    cerr << " input files.\n"
         << "Size of buckets for binning D' values at a particular distance: " << bucket_dist / 1000 << "kb.\n";
    if (dprime_threshold)
        cerr << "D' Threshold set at " << dprime_threshold_level << ". D' values above this limit will be set to 1.0, values below will be set to 0.0.\n";

    //
    // Parse the population map.
    //
    parse_population_map(pmap_path, pop_map, pop_cnts);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<pair<int, string> > files;
    if (!build_file_list(files))
        exit(1);

    cerr << "Identified " << files.size() << " files.\n";

    //
    // Open the log file.
    //
    stringstream log;
    log << "phasedstacks.log";
    string log_path = in_path + log.str();
    ofstream log_fh(log_path.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log_path << "'\n";
        exit(1);
    }

    init_log(log_fh, argc, argv);

    //
    // Load the catalog
    //
    cerr << "Parsing the catalog...\n";
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    bool compressed = false;
    catalog_file << cat_path << "batch_" << batch_id << ".catalog";
    int res = load_loci(catalog_file.str(), catalog, 0, false, compressed);
    if (res == 0) {
        cerr << "Error: Unable to load the catalog '" << catalog_file.str() << "'\n";
        throw exception();
    }
    cerr << "done.\n";

    //
    // Implement the black/white list
    //
    reduce_catalog(catalog, whitelist, blacklist);

    map<int, int> fgt_block_lens, fgt_snp_cnts;
    map<int, int> dp_block_lens, dp_snp_cnts;

    //
    // Vectors to store D' measures of SNPs at bucketed distances.
    //
    vector<double> dprime_buckets, dprime_bucket_cnts;

    for (uint i = 0; i < files.size(); i++) {

        // if (files[i].second != "batch_1.groupV.phase") continue;

        PhasedSummary *psum = NULL;

        if (in_file_type == FileT::fastphase) {
            if ((psum = parse_fastphase(in_path + files[i].second)) == NULL) {
                cerr << "Unable to parse fastPhase input files.\n";
                exit(1);
            }
        } else if (in_file_type == FileT::beagle && haplotypes) {
            if ((psum = parse_beagle_haplotypes(catalog, in_path + files[i].second)) == NULL) {
                cerr << "Unable to parse Beagle input files.\n";
                exit(1);
            }
        } else if (in_file_type == FileT::beagle) {
            if ((psum = parse_beagle(catalog, in_path + files[i].second)) == NULL) {
                cerr << "Unable to parse Beagle input files.\n";
                exit(1);
            }
        }

        //
        // Summarize the genotypes in the populations.
        //
        summarize_phased_genotypes(psum);

        // for (uint j = 0; j < psum->size; j++) {
        //     cerr << "BP: " << psum->nucs[j].bp << "\t"
        //       << "A: "  << std::setw(3) << psum->nucs[j].nuc[0] << " "
        //       << "C: "  << std::setw(3) << psum->nucs[j].nuc[1] << " "
        //       << "G: "  << std::setw(3) << psum->nucs[j].nuc[2] << " "
        //       << "T: "  << std::setw(3) << psum->nucs[j].nuc[3] << "\n";
        // }

        //
        // Calculate D'
        //
        cerr << "Calculating D'...";
        calc_dprime(psum);
        cerr << "done.\n";

        write_dprime(in_path + files[i].second, psum);

        //
        // Generate haplotype blocks based on D'.
        //
        dprime_blocks(in_path + files[i].second, pop_map, psum, dp_block_lens, dp_snp_cnts);

        //
        // Generate haplotype blocks using the four gamete test.
        //
        four_gamete_test(in_path + files[i].second, pop_map, psum, fgt_block_lens, fgt_snp_cnts);

        //
        // Bucket the D' measures by distance between SNPs.
        //
        bucket_dprime(dprime_buckets, dprime_bucket_cnts, psum);

        //
        // Free the Samples objects
        //
        delete psum;
    }

    //
    // Write average D' values bucketed according to their distance in the genome.
    //
    write_buckets(in_path, dprime_buckets, dprime_bucket_cnts);

    //
    // Write the FGT bucketed distances.
    //
    log_fh << "# Distribution of FGT haplotype block lengths.\n";
    map<int, int>::iterator buck_it;
    for (buck_it = fgt_block_lens.begin(); buck_it != fgt_block_lens.end(); buck_it++)
        log_fh << buck_it->first << "\t" << buck_it->second << "\n";

    //
    // Write the FGT bucketed SNP counts.
    //
    log_fh << "\n\n"
           << "# Distribution of FGT SNP counts per haplotype block.\n";
    for (buck_it = fgt_snp_cnts.begin(); buck_it != fgt_snp_cnts.end(); buck_it++)
        log_fh << buck_it->first << "\t" << buck_it->second << "\n";

    //
    // Write the D' haplotype block bucketed distances.
    //
    log_fh << "\n\n"
           << "# Distribution of D' haplotype block lengths.\n";
    for (buck_it = dp_block_lens.begin(); buck_it != dp_block_lens.end(); buck_it++)
        log_fh << buck_it->first << "\t" << buck_it->second << "\n";

    //
    // Write the D' bucketed SNP counts.
    //
    log_fh << "\n\n"
           << "# Distribution of D' SNP counts per haplotype block.\n";
    for (buck_it = dp_snp_cnts.begin(); buck_it != dp_snp_cnts.end(); buck_it++)
        log_fh << buck_it->first << "\t" << buck_it->second << "\n";

    log_fh.close();

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
bucket_dprime(vector<double> &dprime_buckets, vector<double> &dprime_bucket_cnts, PhasedSummary *psum)
{
    uint bucket, dist, max_bucket;
    uint max_dist = 0;

    //
    // Check that we have enough buckets in our vectors. Find the maximum distance between
    // SNPs on this chromosome and add buckets as necessary.
    //
    for (uint i = 0; i < psum->size; i++) {
        for (uint j = i+1; j < psum->size; j++) {

            if (psum->nucs[i].freq < minor_freq_lim ||
                psum->nucs[j].freq < minor_freq_lim)
                continue;

            if (write_zeros == false && psum->dprime[i][j].chisq_p == false)
                continue;

            dist     = psum->nucs[j].bp - psum->nucs[i].bp;
            max_dist = dist > max_dist ? dist : max_dist;
        }
    }
    max_bucket = max_dist / bucket_dist;
    if (dprime_buckets.size() < max_bucket) {
        uint cnt = max_bucket + 1 - dprime_buckets.size();
        for (uint i = 0; i < cnt; i++) {
            dprime_buckets.push_back(0.0);
            dprime_bucket_cnts.push_back(0.0);
        }
    }

    //
    // Populate buckets
    //
    for (uint i = 0; i < psum->size; i++) {
        for (uint j = i+1; j < psum->size; j++) {

            if (psum->nucs[i].freq < minor_freq_lim ||
                psum->nucs[j].freq < minor_freq_lim)
                continue;

            if (write_zeros == false && psum->dprime[i][j].chisq_p == false)
                continue;

            bucket = ((psum->nucs[j].bp - psum->nucs[i].bp) / bucket_dist);

            dprime_buckets[bucket] += (psum->dprime[i][j].chisq_p ? psum->dprime[i][j].dprime : 0.0);
            dprime_bucket_cnts[bucket]++;
        }
    }

    return 0;
}

int
write_buckets(string path, vector<double> &dprime_buckets, vector<double> &dprime_bucket_cnts)
{
    //
    // Write the bucketed D' data for plotting.
    //
    stringstream file;
    file << path << "Dprime_dist_buckets" << bucket_dist/1000 << "kb.tsv";

    cerr << "Writing bucketed D' data to '" << file.str() << "'...";

    ofstream fh(file.str().c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening D' file '" << file.str() << "'\n";
        exit(1);
    }

    fh << "# Distance (Kb)\tD' Average\n";

    for (uint i = 0; i < dprime_buckets.size(); i++)
        fh << (i * bucket_dist) << "\t"
           << std::setprecision(3) << (dprime_buckets[i] / dprime_bucket_cnts[i]) << "\n";

    fh.close();

    cerr << "done\n";

    return 0;
}

int
four_gamete_test(string path, map<string, int> &pop_map, PhasedSummary *psum, map<int, int> &len_buckets, map<int, int> &snp_buckets)
{
    //
    // Write haplotypes as found by the four gamete test:
    //   Wang, et al., Am. J. Hum. Genet. 71:1227–1234, 2002
    //
    string file = path + ".fgt.tsv";

    cerr << "Determining four gamete test haplotypes blocks, writing to:\n    '" << file << "'...\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening FGT file '" << file << "'\n";
        exit(1);
    }

    fh << "# ID\tStart\tEnd\tLen\tSNP Count\tHaplotype Count\tHaplotype\tPopulations\tHapPopCnt\n";

    uint id = 1;
    uint start, end=-1, cnt, dist;
    bool bound;
    map<int, int> buckets, snps;

    for (uint i = 0; i < psum->size; i++) {
        if (psum->nucs[i].freq < minor_freq_lim)
            continue;

        //
        // Start a new block.
        //
        start  = i;
        bound  = false;
        cnt    = 0;
        uint j = i;

        do {
            if (psum->nucs[j].freq < minor_freq_lim) {
                j++;
                continue;
            }

            for (int k = j; k >= (int) start; k--) {

                if (psum->nucs[k].freq < minor_freq_lim)
                    continue;

                if (psum->recomb[k][j] == true) {
                    bound = true;
                    end   = j;
                }
            }

            j++;
            cnt++;
        } while (bound == false && j < psum->size);

        if (j == psum->size)
            end = j - 1;

        fh << id << "\t"
           << psum->nucs[start].bp << "\t"
           << psum->nucs[end].bp << "\t"
           << psum->nucs[end].bp - psum->nucs[start].bp + 1 << "\t"
           << cnt << "\t";
        //
        // Bucket the SNP counts for plotting.
        //
        snps[cnt]++;

        //
        // Bucket the haplotype block lengths for plotting.
        //
        dist = (psum->nucs[end].bp - psum->nucs[start].bp + 1) / 10000 * 10000;
        buckets[dist]++;

        enumerate_haplotypes(fh, pop_map, psum, start, end);

        i = end;
        id++;
    }

    //
    // Write the bucketed distances.
    //
    fh << "\n\n"
       << "# Distribution of FGT haplotype block lengths.\n";
    map<int, int>::iterator it;
    for (it = buckets.begin(); it != buckets.end(); it++) {
        fh << it->first << "\t" << it->second << "\n";

        len_buckets[it->first] += it->second;
    }

    //
    // Write the bucketed SNP counts.
    //
    fh << "\n\n"
       << "# Distribution of SNP counts per FGT haplotype block.\n";
    for (it = snps.begin(); it != snps.end(); it++) {
        fh << it->first << "\t" << it->second << "\n";

        snp_buckets[it->first] += it->second;
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
dprime_blocks(string path, map<string, int> &pop_map, PhasedSummary *psum, map<int, int> &len_buckets, map<int, int> &snp_buckets)
{
    //
    // Generate haplotype blocks according to strength of linkage disequilibrium measured using D'.
    //   Stacey B. Gabriel et al. (2002). The Structure of Haplotype Blocks in the Human Genome. Science 296:2225-2229
    //
    string file = path + ".dpblocks.tsv";

    cerr << "Determining D' haplotypes blocks, writing to:\n    '" << file << "'...\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening D' blocks file '" << file << "'\n";
        exit(1);
    }

    fh << "# ID\tStart\tEnd\tLen\tSNP Count\tHaplotype Count\tHaplotype\tPopulations\tHapPopCnt\n";

    uint dist;
    set<int> loci;
    vector<pair<int, int> > ld_pairs;
    map<int, vector<int> > ld_map;
    map<int, int> buckets, snps;

    uint tot_pairs    = 0;
    uint recomb_pairs = 0;

    for (uint i = 0; i < psum->size; i++) {
        if (psum->nucs[i].freq < minor_freq_lim)
            continue;

        for (uint j = i+1; j < psum->size; j++) {
            if (psum->nucs[j].freq < minor_freq_lim)
                continue;

            tot_pairs++;
            dist = psum->nucs[j].bp - psum->nucs[i].bp + 1;

            //
            // Does this pair of markers show a strong measure of LD?
            //
            if (psum->dprime[i][j].ci_high > 0.98 &&
                psum->dprime[i][j].ci_low  > 0.7 &&
                dist <= max_pair_dist) {
                psum->dprime[i][j].type = strong_ld;

                ld_pairs.push_back(make_pair(i, j));
                ld_map[i].push_back(j);

                loci.insert(i);
                loci.insert(j);
            }

            //
            // Does this pair of markers show a strong measure of historical recombination?
            //
            if (psum->dprime[i][j].ci_high < 0.9) {
                psum->dprime[i][j].type = recomb;
                recomb_pairs++;
            }
        }
    }

    // map<int, vector<int> >::iterator it;
    // for (it = ld_map.begin(); it != ld_map.end(); it++) {
    //  cerr << "      " << it->first << " ->\n";
    //  for (uint i = 0; i < it->second.size(); i++)
    //      cerr << "         " << it->second[i] << " dist: " << (psum->nucs[it->second[i]].bp - psum->nucs[it->first].bp + 1) << "bp\n";
    // }

    cerr << "    Total pairs examined: " << tot_pairs
         << "; Strong LD pairs: " << ld_pairs.size()
         << "; Recombination pairs: " << recomb_pairs
         << "; Informative markers: " << std::setprecision(3)
         << ((double) (ld_pairs.size() + recomb_pairs) / (double) tot_pairs) * 100 << "%\n";

    //
    // Convert our list of loci into an ordered, linked list, where each node
    // represents a haplotype block.
    //
    dPrimeBlocks blocks;
    blocks.initialize(loci);

    //
    // Merge nodes together where D' is strong enough to maintain the block.
    //
    HBlock *cur;

    cur = blocks.head();

    do {
        //
        // Can we merge these two nodes together?
        //
        if (check_adjacent_blocks(psum, cur)) {
            // cerr << "    Merging blocks: ";
            // for (uint i = 0; i < cur->loci.size(); i++)
            //  cerr << cur->loci[i] << ", ";
            // cerr << " and ";
            // for (uint i = 0; i < cur->next->loci.size(); i++)
            //  cerr << cur->next->loci[i] << ", ";
            // cerr << "\n";
            blocks.merge_adjacent(cur);
        } else {
            cur = cur->next;
        }
    } while (cur->next != NULL);

    // blocks.print();

    //
    // Write the blocks.
    //
    uint id = 1;
    uint start, end;

    cur = blocks.head();

    do {
        start = *(cur->loci.begin());
        end   = *(cur->loci.rbegin());

        fh << id << "\t"
           << psum->nucs[start].bp << "\t"
           << psum->nucs[end].bp << "\t"
           << psum->nucs[end].bp - psum->nucs[start].bp + 1 << "\t"
           << cur->loci.size() << "\t";

        //
        // Bucket the SNP counts for plotting.
        //
        snps[cur->loci.size()]++;

        //
        // Bucket the haplotype block lengths for plotting.
        //
        dist = (psum->nucs[end].bp - psum->nucs[start].bp + 1) / 10000 * 10000;
        buckets[dist]++;

        enumerate_haplotypes(fh, pop_map, psum, start, end);

        id++;
        cur = cur->next;
    } while (cur != NULL);

    //
    // Write the bucketed distances.
    //
    fh << "\n\n"
       << "# Distribution of D' haplotype block lengths.\n";
    map<int, int>::iterator it;
    for (it = buckets.begin(); it != buckets.end(); it++) {
        fh << it->first << "\t" << it->second << "\n";

        len_buckets[it->first] += it->second;
    }

    //
    // Write the bucketed SNP counts.
    //
    fh << "\n\n"
       << "# Distribution of SNP counts per D' haplotype block.\n";
    for (it = snps.begin(); it != snps.end(); it++) {
        fh << it->first << "\t" << it->second << "\n";

        snp_buckets[it->first] += it->second;
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

bool
check_adjacent_blocks(PhasedSummary *psum, HBlock *block)
{
    //
    // Create a list of all loci contained in the two blocks.
    //
    uint start = *(block->loci.begin());
    uint end   = *(block->next->loci.rbegin());

    //
    // Check the D' measure between each pair in the proposed combined block.
    //
    double tot       = 0.0;
    double strong_ld = 0.0;

    for (uint i = start; i <= end; i++) {
        if (psum->nucs[i].freq < minor_freq_lim)
            continue;

        for (uint j = i + 1; j <= end; j++) {
            if (psum->dprime[i][j].type == uninformative ||
                psum->nucs[j].freq < minor_freq_lim)
                continue;

            tot++;
            if (psum->dprime[i][j].type == strong_ld)
                strong_ld++;
        }
    }

    // cerr << "Comparing range " << start << " to " << end
    //   << "; total pairs: " << tot << "; strong LD: " << strong_ld
    //   << "; proportion: " << std::setprecision(3) << strong_ld / tot << "\n";

    if (strong_ld / tot >= min_inform_pairs)
        return true;

    return false;
}

HBlock *
dPrimeBlocks::merge_adjacent(HBlock *a)
{
    //
    // Merge two adjacent nodes.
    //
    HBlock *b = a->next;

    for (uint i = 0; i < b->loci.size(); i++)
        a->loci.push_back(b->loci[i]);
    a->next = b->next;
    delete b;
    return a;
}

HBlock *
dPrimeBlocks::initialize(set<int> &loci)
{
    set<int>::iterator it, prev_it;
    HBlock *cur, *next;

    this->_head = new HBlock;
    it = loci.begin();
    this->_head->loci.push_back(*it);
    it++;

    // //
    // // Create a node from each locus and add to it all immediately adjacent loci.
    // //
    // do {
    //  this->_head->loci.push_back(*it);
    //  prev_it = it;
    //  it++;
    // } while (it != loci.end() && (*prev_it) + 1 == *it);

    next = this->_head;

    // if (it == loci.end()) return this->_head;

    do {
        cur = new HBlock;
        cur->loci.push_back(*it);
        it++;

        // do {
        //     cur->loci.push_back(*it);
        //     prev_it = it;
        //     it++;
        // } while (it != loci.end() &&
        //       (*prev_it) + 1 == *it);

        next->next = cur;
        next = next->next;

    } while (it != loci.end());

    return this->_head;
}

int
dPrimeBlocks::print()
{
    HBlock *cur = this->_head;
    while (cur != NULL) {
        for (uint i = 0; i < cur->loci.size(); i++) {
            if (i > 0)
                cerr << ", ";
            cerr << cur->loci[i];
        }
        cerr << "\n";

        cur = cur->next;
    }

    return 0;
}

int
enumerate_haplotypes(ofstream &fh, map<string, int> &pop_map, PhasedSummary *psum, uint start, uint end)
{
    map<string, map<int, int> >::iterator it;
    map<string, map<int, int> > haplotypes;
    map<int, int>::iterator sit;
    string haplotype;
    set<int> pops;

    //
    // Enumerate all haplotypes occurring in this block.
    //
    for (uint k = 0; k < psum->sample_cnt; k++) {

        for (uint n = start; n <= end; n++)
            if (psum->nucs[n].freq >= minor_freq_lim)
                haplotype += psum->samples[k].nucs_1[n];

        pops.insert(pop_map[psum->samples[k].name]);

        if (haplotypes.count(haplotype) == 0)
            haplotypes[haplotype][pop_map[psum->samples[k].name]] = 1;
        else
            haplotypes[haplotype][pop_map[psum->samples[k].name]]++;
        haplotype.clear();
    }

    for (uint k = 0; k < psum->sample_cnt; k++) {

        for (uint n = start; n <= end; n++)
            if (psum->nucs[n].freq >= minor_freq_lim)
                haplotype += psum->samples[k].nucs_2[n];

        pops.insert(pop_map[psum->samples[k].name]);

        if (haplotypes.count(haplotype) == 0)
            haplotypes[haplotype][pop_map[psum->samples[k].name]] = 1;
        else
            haplotypes[haplotype][pop_map[psum->samples[k].name]]++;
        haplotype.clear();
    }

    //
    // Write the haplotypes.
    //
    float tot = 0.0;

    fh << haplotypes.size() << "\t";
    for (it = haplotypes.begin(); it != haplotypes.end(); it++) {
        //
        // Haplotypes are stored per population; sum them up here.
        //
        for (sit = it->second.begin(); sit != it->second.end(); sit++)
            tot += sit->second;

        if (it != haplotypes.begin())
            fh << ",";
        fh << it->first << "|"
           << std::setprecision(3) << tot / ((float) psum->sample_cnt * 2.0);
    }
    fh << "\t";

    set<int>::iterator pit;
    stringstream       pops_str;

    //
    // Write which populations this haplotype block occurs in.
    //
    if (pops.size() == 0)
        fh << "-1\t";
    else
        for (pit = pops.begin(); pit != pops.end(); pit++)
            pops_str << *pit << ",";
    fh << pops_str.str().substr(0, pops_str.str().length()-1);
    pops_str.str("");

    //
    // Write the frequency of occurence of each haplotype in each population.
    //
    for (it = haplotypes.begin(); it != haplotypes.end(); it++) {
        pops_str << "\t";
        for (pit = pops.begin(); pit != pops.end(); pit++)
            pops_str << (it->second)[*pit] << "|"
                     << std::setprecision(3)
                     << (float) (it->second)[*pit] / (float) (pop_cnts[*pit] * 2.0) << ",";
        fh << pops_str.str().substr(0, pops_str.str().length()-1);
        pops_str.str("");
    }
    fh << "\n";

    return 0;
}

int
calc_dprime(PhasedSummary *psum)
{
    #pragma omp parallel
    {
        char   allele_A, allele_a, allele_B, allele_b;
        double freq_A,     freq_a,   freq_B,   freq_b;
        double freq_AB,   freq_Ab,  freq_aB,  freq_ab;
        double D, min, var, chisq;
        double tot = psum->sample_cnt  * 2.0;
        uint   hap_cnt;

        #pragma omp for schedule(dynamic, 1)
        for (uint i = 0; i < psum->size; i++) {
            //
            // Assign nucleotides to allele A, and a.
            //
            assign_alleles(psum->nucs[i], allele_A, allele_a, freq_A, freq_a);

            for (uint j = i+1; j < psum->size; j++) {
                //
                // Assign nucleotides to allele B, and b.
                //
                assign_alleles(psum->nucs[j], allele_B, allele_b, freq_B, freq_b);

                freq_AB = 0.0;
                freq_Ab = 0.0;
                freq_aB = 0.0;
                freq_ab = 0.0;
                hap_cnt = 0;
                D       = 0.0;

                //
                // Tally up haplotype frequencies.
                //
                for (uint k = 0; k < psum->sample_cnt; k++) {

                    if (psum->samples[k].nucs_1[i] == allele_A &&
                        psum->samples[k].nucs_1[j] == allele_B)
                        freq_AB++;
                    else if (psum->samples[k].nucs_1[i] == allele_A &&
                             psum->samples[k].nucs_1[j] == allele_b)
                        freq_Ab++;
                    else if (psum->samples[k].nucs_1[i] == allele_a &&
                             psum->samples[k].nucs_1[j] == allele_B)
                        freq_aB++;
                    else if (psum->samples[k].nucs_1[i] == allele_a &&
                             psum->samples[k].nucs_1[j] == allele_b)
                        freq_ab++;

                    if (psum->samples[k].nucs_2[i] == allele_A &&
                        psum->samples[k].nucs_2[j] == allele_B)
                        freq_AB++;
                    else if (psum->samples[k].nucs_2[i] == allele_A &&
                             psum->samples[k].nucs_2[j] == allele_b)
                        freq_Ab++;
                    else if (psum->samples[k].nucs_2[i] == allele_a &&
                             psum->samples[k].nucs_2[j] == allele_B)
                        freq_aB++;
                    else if (psum->samples[k].nucs_2[i] == allele_a &&
                             psum->samples[k].nucs_2[j] == allele_b)
                        freq_ab++;
                }

                freq_AB = freq_AB / tot;
                freq_Ab = freq_Ab / tot;
                freq_aB = freq_aB / tot;
                freq_ab = freq_ab / tot;

                //
                // Using the four-gamete test, check whether recombination has occurred
                // between these two loci.
                //  Four-gamete test: if no recombination has occurred between any two loci (SNPs) there will
                //  be three haplotypes present, if recombination has occurred there will be four haplotypes.
                //
                hap_cnt += freq_AB > 0 ? 1 : 0;
                hap_cnt += freq_Ab > 0 ? 1 : 0;
                hap_cnt += freq_aB > 0 ? 1 : 0;
                hap_cnt += freq_ab > 0 ? 1 : 0;

                if (hap_cnt == 3)
                    psum->recomb[i][j] = false;
                else
                    psum->recomb[i][j] = true;


                D = freq_AB - (freq_A * freq_B);
                // cerr << "D_AB: " << D << "; ";
                // D = freq_Ab - (freq_A * freq_b);
                // cerr << "D_Ab: " << D << "; ";
                // D = freq_aB - (freq_a * freq_B);
                // cerr << "D_aB: " << D << "; ";
                // D = freq_ab - (freq_a * freq_b);
                // cerr << "D_ab: " << D << "\n";
                // cerr << "    freq_AB: " << freq_AB << "; freq_Ab: " << freq_Ab << "; freq_aB: " << freq_aB << "; freq_ab: " << freq_ab << "\n";

                if (D > 0) {
                    min = (freq_A * freq_b) < (freq_a * freq_B) ? (freq_A * freq_b) : (freq_a * freq_B);
                    psum->dprime[i][j].dprime = min == 0 ? 0.0 : D / min;
                } else {
                    min = (freq_A * freq_B) < (freq_a * freq_b) ? (freq_A * freq_B) : (freq_a * freq_b);
                    psum->dprime[i][j].dprime = min == 0 ? 0.0 :(-1 * D) / min;
                }

                //
                // Test D against a chi square distribution with 1 degree of freedom to show
                // whether these two loci have a D that is statistically significantly different from 0.
                //
                chisq = (tot * (D * D)) / (freq_A * freq_a * freq_B * freq_b);
                if (chisq >= chi_sq_limit)
                    psum->dprime[i][j].chisq_p = true;

                //
                // Calculate variance and confidence limits.
                //
                if (psum->dprime[i][j].chisq_p) {
                    var   = (1.0 / tot) * ((freq_A * freq_a * freq_B * freq_b) + ((1 - (2 * freq_A)) * (1 - (2 * freq_B)) * D) - (D * D));
                    psum->dprime[i][j].var = var;
                    psum->dprime[i][j].ci_high = psum->dprime[i][j].dprime + (1.96 * sqrt(var));
                    psum->dprime[i][j].ci_low  = psum->dprime[i][j].dprime - (1.96 * sqrt(var));
                } else {
                    psum->dprime[i][j].ci_high = 0.0;
                    psum->dprime[i][j].ci_low  = 0.0;
                }
            }
        }
    }

    return 0;
}

int
assign_alleles(NucSum nsum, char &p_allele, char &q_allele, double &p_freq, double &q_freq)
{
    p_allele = 0;
    q_allele = 0;

    uint  i   = 0;
    float tot = 0;

    while (p_allele == 0 && i < 4) {
        if (nsum.nuc[i] > 0) {
            tot += nsum.nuc[i];
            switch(i) {
            case 0:
                p_allele = 'A';
                p_freq   = nsum.nuc[0];
                break;
            case 1:
                p_allele = 'C';
                p_freq   = nsum.nuc[1];
                break;
            case 2:
                p_allele = 'G';
                p_freq   = nsum.nuc[2];
                break;
            case 3:
                p_allele = 'T';
                p_freq   = nsum.nuc[3];
                break;
            }
        }
        i++;
    }
    while (q_allele == 0 && i < 4) {
        if (nsum.nuc[i] > 0) {
            tot += nsum.nuc[i];
            switch(i) {
            case 1:
                q_allele = 'C';
                q_freq   = nsum.nuc[1];
                break;
            case 2:
                q_allele = 'G';
                q_freq   = nsum.nuc[2];
                break;
            case 3:
                q_allele = 'T';
                q_freq   = nsum.nuc[3];
                break;
            }
        }
        i++;
    }

    p_freq = p_freq / tot;
    q_freq = 1 - p_freq;

    return 0;
}

int
write_dprime(string path, PhasedSummary *psum)
{
    //
    // Write the D' data for plotting as a heatmap.
    //
    string file = path + ".dprime.tsv";

    cerr << "Writing D' data to '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening D' file '" << file << "'\n";
        exit(1);
    }

    fh << "# Basepair 1\tBasepair 2\tD'\tCorrected D'\tVariance\tCI Low\tCI High\n";

    double dprime = 0.0;

    for (uint i = 0; i < psum->size; i++) {
        for (uint j = i+1; j < psum->size; j++) {

            if (psum->nucs[i].freq < minor_freq_lim ||
                psum->nucs[j].freq < minor_freq_lim)
                continue;

            dprime = psum->dprime[i][j].dprime;

            if (dprime_threshold)
                dprime = dprime >= dprime_threshold_level ? 1.0 : 0.0;

            if (write_zeros == false && (dprime == 0.0 || psum->dprime[i][j].chisq_p == false))
                continue;

            fh << psum->nucs[i].bp << "\t"
               << psum->nucs[j].bp << "\t"
               << std::setprecision(3) << dprime << "\t"
               << std::setprecision(3) << (psum->dprime[i][j].chisq_p ? dprime : 0.0) << "\t"
               << psum->dprime[i][j].var << "\t"
               << psum->dprime[i][j].ci_low  << "\t"
               << psum->dprime[i][j].ci_high << "\n";
        }
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
summarize_phased_genotypes(PhasedSummary *psum)
{
    //
    // Construct a two dimensional array out of all the nucleotide arrays in the samples.
    //
    char **gtypes = new char *[psum->sample_cnt];

    for (uint i = 0; i < psum->sample_cnt; i++) {
        gtypes[i] = psum->samples[i].nucs_1;
    }

    //
    // Sum up the occurences of each nucleotide.
    //
    for (uint i = 0; i < psum->size; i++) {
        for (uint j = 0; j < psum->sample_cnt; j++) {
            switch(gtypes[j][i]) {
            case 'A':
                psum->nucs[i].nuc[0]++;
                break;
            case 'C':
                psum->nucs[i].nuc[1]++;
                break;
            case 'G':
                psum->nucs[i].nuc[2]++;
                break;
            case 'T':
                psum->nucs[i].nuc[3]++;
                break;
            case 'N':
            default:
                break;
            }
        }
    }

    //
    // Repeat for the second set of phased genotypes.
    //
    for (uint i = 0; i < psum->sample_cnt; i++) {
        gtypes[i] = psum->samples[i].nucs_2;
    }

    //
    // Sum up the occurences of each nucleotide.
    //
    for (uint i = 0; i < psum->size; i++) {
        for (uint j = 0; j < psum->sample_cnt; j++) {
            switch(gtypes[j][i]) {
            case 'A':
                psum->nucs[i].nuc[0]++;
                break;
            case 'C':
                psum->nucs[i].nuc[1]++;
                break;
            case 'G':
                psum->nucs[i].nuc[2]++;
                break;
            case 'T':
                psum->nucs[i].nuc[3]++;
                break;
            case 'N':
            default:
                break;
            }
        }

        //
        // Calculate minor allele frequency.
        //
        float tot  = (float) psum->sample_cnt * 2.0;
        float freq = 0.0;
        for (uint j = 0; j < 4; j++) {
            if (psum->nucs[i].nuc[j] > 0) {
                freq = (float) psum->nucs[i].nuc[j] / tot;
                psum->nucs[i].freq = freq < psum->nucs[i].freq ? freq : psum->nucs[i].freq;
            }
        }
    }

    delete [] gtypes;

    return 0;
}

//
// Code to parse fastPhase format.
//
PhasedSummary *
parse_fastphase(string path)
{
    ifstream    fh;
    char        line[max_len];
    string      buf, filepath;
    const char *p, *q, *end;
    int         i, sindex, pos;

    memset(line, '\0', max_len);

    //
    // Read in the original fastPhase export from Stacks to obtain the original base pair positions.
    //
    //
    // Open the file for reading
    //
    filepath = path + ".inp";
    fh.open(filepath.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening input file '" << path << "'\n";
        return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    int  num_samples, num_genotypes;
    char bp[id_len];

    //
    // Get the number of samples in the dataset.
    //
    fh.getline(line, max_len);
    num_samples = is_integer(line);

    if (num_samples < 0) {
        cerr << "Unable to find the number of samples, should be the first line.\n";
        return NULL;
    }

    //
    // Get the number of genotypes in the dataset.
    //
    fh.getline(line, max_len);
    num_genotypes = is_integer(line);

    if (num_genotypes < 0) {
        cerr << "Unable to find the number of genotypes, should be the second line.\n";
        return NULL;
    }

    PhasedSummary *psum = new PhasedSummary(num_samples, num_genotypes);

    //
    // Get the set of base pair positions.
    //
    buf.clear();
    do {
        fh.clear();
        fh.getline(line, max_len);
        buf += line;
    } while (fh.fail() && !fh.bad() && !fh.eof());

    i   = 0;
    p   = buf.c_str();
    end = p + buf.length();

    if (*p != 'P') {
        cerr << "Unable to locate line of basepair positions, should be the third line.\n";
        delete psum;
        return NULL;
    }
    for (p += 2, q = p; p < end; p++, q++) {
        while (*q != ' ' && q < end) {
            q++;
        }
        strncpy(bp, p, q - p);
        bp[q - p] = '\0';
        pos = is_integer(bp);

        if (pos < 0) {
            cerr << "Unable to parse base pair positions.\n";
            delete psum;
            return NULL;
        } else {
            psum->nucs[i].bp = (uint) pos;
        }

        i++;
        p = q;
    }

    fh.close();

    //
    // Open the file for reading
    //
    filepath = path + "_hapguess_switch.out";
    fh.open(filepath.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening input file '" << path << "'\n";
        return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    //
    // Read from the "*_hapguess_switch.out" file until we hit the genotypes section
    // marked by the string "BEGIN GENOTYPES".
    //
    do {
        fh.getline(line, max_len);

        if (!fh.good()) {
            cerr << "Unable to find file section entitled 'BEGIN GENOTYPES'\n";
            delete psum;
            return NULL;
        }

    } while (strcmp(line, "BEGIN GENOTYPES") != 0);

    //
    // Now read lines from the file in groups of three:
    //   1. Sample label
    //   2. Phased genotypes from chromosome 1
    //   3. Phased genotypes from chromosome 2
    // Stop reading individuals when we encounter the string, "END GENOTYPES".
    //
    fh.getline(line, max_len);

    do {
        //
        // Create a new Sample object and store the sample label.
        //
        sindex = psum->add_sample(line);

        //
        // Get the first set of phased genotypes.
        //
        buf.clear();
        do {
            fh.clear();
            fh.getline(line, max_len);
            buf += line;
        } while (fh.fail() && !fh.bad() && !fh.eof());

        //
        // Count the number of genotypes on this line (they should be space deliniated).
        //
        i = 0;
        for (p = buf.c_str(); *p != '\0'; p++)
            if (*p != ' ') psum->samples[sindex].size++;
        //
        // Store the genotypes into our internal buffer.
        //
        psum->samples[sindex].nucs_1 = new char[psum->samples[sindex].size];
        for (p = buf.c_str(); *p != '\0'; p++) {
            if (*p == ' ') continue;
            psum->samples[sindex].nucs_1[i] = *p;
            i++;
        }

        // len = strlen(line);
        // if (line[len - 1] == '\r') line[len - 1] = '\0';

        //
        // Get the second set of phased genotypes.
        //
        buf.clear();
        do {
            fh.clear();
            fh.getline(line, max_len);
            buf += line;
        } while (fh.fail() && !fh.bad() && !fh.eof());

        i = 0;
        psum->samples[sindex].nucs_2 = new char[psum->samples[sindex].size];
        for (p = buf.c_str(); *p != '\0'; p++) {
            if (*p == ' ') continue;
            psum->samples[sindex].nucs_2[i] = *p;
            i++;
        }

        //
        // Get the sample label of the next record.
        //
        fh.getline(line, max_len);

    } while (strcmp(line, "END GENOTYPES") != 0 && fh.good());

    fh.close();

    return psum;
}

//
// Code to parse Beagle format.
//
PhasedSummary *
parse_beagle(map<int, CSLocus *> &catalog, string path)
{
    gzFile      gz_fh;
    char        *line;
    string      buf, filepath;
    const char *p, *q;
    uint        len, line_len, i, sindex;
    bool        eol;

    line_len = max_len;
    line = new char[line_len];
    memset(line, '\0', line_len);

    //
    // Open the Beagle file for reading
    //
    filepath = path + ".phased.gz";
    gz_fh = gzopen(filepath.c_str(), "rb");
    if (!gz_fh) {
        cerr << "Failed to open gzipped file '" << filepath << "': " << strerror(errno) << ".\n";
        return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    vector<string> parts;
    uint num_samples   = 0;
    uint num_genotypes = 0;
    char cat_loc_str[id_len], col_str[id_len];

    //
    // Parse the file twice. On the first round:
    //   1. Determine the number of samples in the dataset (column count)
    //   2. Determine the number of markers (row count).
    // On the second round, parse the SNP genotypes.
    //

    //
    // Read each line in the file. If it starts with:
    //   '#' it is a comment, skip the line.
    //   'I' it is the list of samples, parse them.
    //   'S' is the population ID for each SNP, skip this line.
    //   'M' is a marker, count the number of markers.
    //
    do {
        eol = false;
        buf.clear();
        do {
            gzgets(gz_fh, line, line_len);
            buf += line;

            len = strlen(line);
            if (len > 0 && line[len - 1] == '\n') {
                eol = true;
                line[len - 1] = '\0';
            }
        } while (!gzeof(gz_fh) && !eol);

        if (line_len < buf.length()) {
            // cerr << "Resizing line buffer from " << line_len << " to " << buf.length() << "\n";
            delete [] line;
            line = new char[buf.length() + 1];
            line_len = buf.length() + 1;
            memset(line, '\0', line_len);
        }

        if (buf[0] == 'M') {
            num_genotypes++;
        } else if (buf[0] == 'I') {
            //
            // Count the number of samples.
            //
            parse_ssv(buf.c_str(), parts);
            num_samples = (parts.size() - 2) / 2;
        }

    } while (!gzeof(gz_fh));

    PhasedSummary *psum = new PhasedSummary(num_samples, num_genotypes);

    for (uint j = 2; j < parts.size(); j++) {
        if (j % 2 == 0) {
            sindex = psum->add_sample(parts[j]);
            psum->samples[sindex].size   = num_genotypes;
            psum->samples[sindex].nucs_1 = new char[psum->samples[sindex].size];
            psum->samples[sindex].nucs_2 = new char[psum->samples[sindex].size];
        }
    }

    cerr << "  Found " << num_samples << " samples; " << num_genotypes << " genotypes.\n";

    gzrewind(gz_fh);

    uint marker_num = 0;
    memset(line, '\0', line_len);

    do {
        do {
            gzgets(gz_fh, line, line_len);
        } while (!gzeof(gz_fh) && line[0] != 'M');

        len = strlen(line);

        if (len == 0) break;
        if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';

        parse_ssv(line, parts);

        //
        // Parse the catalog locus ID and the column number of the SNP:
        //   e.g. LocId_column or 10329_37
        //
        p = parts[1].c_str();
        for (q = p + 1; *q != '_' && *q != '\0'; q++);
        strncpy(cat_loc_str, p, q - p);
        cat_loc_str[q-p] = '\0';
        q++;
        strcpy(col_str, q);

        psum->nucs[marker_num].clocus = is_integer(cat_loc_str);
        psum->nucs[marker_num].col    = is_integer(col_str);

        //
        // Store the genotypes into our internal buffer.
        //
        sindex = 0;
        i      = 2;
        while (i < parts.size()) {
            p = parts[i].c_str();
            psum->samples[sindex].nucs_1[marker_num] = *p;
            i++;
            p = parts[i].c_str();
            psum->samples[sindex].nucs_2[marker_num] = *p;
            i++;
            sindex++;
        }

        marker_num++;

    } while (!gzeof(gz_fh));

    gzclose(gz_fh);

    //
    // Use the catalog to look up the basepair positions for each catalog locus.
    //
    CSLocus *loc;
    for (i = 0; i < psum->size; i++) {
        loc = catalog[psum->nucs[i].clocus];
        psum->nucs[i].bp = loc->sort_bp(psum->nucs[i].col);
    }

    return psum;
}

//
// Code to parse Beagle format.
//
PhasedSummary *
parse_beagle_haplotypes(map<int, CSLocus *> &catalog, string path)
{
    gzFile      gz_fh;
    char        *line;
    string      buf, filepath;
    const char *p;
    uint        len, line_len, i, j, sindex;
    bool        eol;

    line_len = max_len;
    line = new char[line_len];
    memset(line, '\0', line_len);

    //
    // Open the Beagle file for reading
    //
    filepath = path + ".phased.gz";
    gz_fh = gzopen(filepath.c_str(), "rb");
    if (!gz_fh) {
        cerr << "Failed to open gzipped file '" << filepath << "': " << strerror(errno) << ".\n";
        return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    vector<string> parts, samples;
    uint num_samples   = 0;
    uint num_genotypes = 0;
    uint cat_loc;

    //
    // Parse the file twice. On the first round:
    //   1. Determine the number of samples in the dataset (column count)
    //   2. Determine the number of markers (row count).
    // On the second round, parse the SNP genotypes.
    //

    //
    // Read each line in the file. If it starts with:
    //   '#' it is a comment, skip the line.
    //   'I' it is the list of samples, parse them.
    //   'S' is the population ID for each SNP, skip this line.
    //   'M' is a marker, count the number of markers.
    //
    do {
        eol = false;
        buf.clear();
        do {
            gzgets(gz_fh, line, line_len);
            buf += line;

            len = strlen(line);
            if (len > 0 && line[len - 1] == '\n') {
                eol = true;
                line[len - 1] = '\0';
            }
        } while (!gzeof(gz_fh) && !eol);

        if (line_len < buf.length()) {
            // cerr << "Resizing line buffer from " << line_len << " to " << buf.length() << "\n";
            delete [] line;
            line = new char[buf.length() + 1];
            line_len = buf.length() + 1;
            memset(line, '\0', line_len);
        }

        if (buf[0] == 'M') {
            //
            // Count the number of genotypes by counting the number or nucleotides in each
            // haplotype for each marker.
            //
            parse_ssv(buf.c_str(), parts);
            num_genotypes += parts[2].length();

        } else if (buf[0] == 'I') {
            //
            // Count the number of samples.
            //
            parse_ssv(buf.c_str(), samples);
            num_samples = (samples.size() - 2) / 2;
        }

    } while (!gzeof(gz_fh));

    PhasedSummary *psum = new PhasedSummary(num_samples, num_genotypes);

    for (uint j = 2; j < samples.size(); j++) {
        if (j % 2 == 0) {
            sindex = psum->add_sample(samples[j]);
            psum->samples[sindex].size   = num_genotypes;
            psum->samples[sindex].nucs_1 = new char[psum->samples[sindex].size];
            psum->samples[sindex].nucs_2 = new char[psum->samples[sindex].size];
        }
    }

    cerr << "  Found " << num_samples << " samples; " << num_genotypes << " genotypes.\n";

    gzrewind(gz_fh);

    CSLocus *loc;
    uint hap_len    = 0;
    uint marker_num = 0;
    memset(line, '\0', line_len);

    do {
        do {
            gzgets(gz_fh, line, line_len);
        } while (!gzeof(gz_fh) && line[0] != 'M');

        len = strlen(line);

        if (len == 0) break;
        if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';

        parse_ssv(line, parts);

        //
        // Use the catalog to look up the basepair positions for each catalog locus.
        //
        cat_loc = is_integer(parts[1].c_str());
        loc     = catalog[cat_loc];
        hap_len = parts[2].length();

        if (hap_len != loc->snps.size())
            cerr << "Haplotypes don't match between catalog and beagle; Locus ID: " << loc->id << "; beagle hap len: " << hap_len << "; catalog hap len: " << loc->snps.size() << "\n";

        for (j = 0, i = marker_num; i < marker_num + hap_len; i++, j++) {
            psum->nucs[i].clocus = cat_loc;
            psum->nucs[i].col    = loc->snps[j]->col;
            psum->nucs[i].bp     = loc->sort_bp(psum->nucs[i].col);
        }

        //
        // Store the genotypes into our internal buffer.
        //
        sindex = 0;
        i      = 2;
        while (i < parts.size()) {
            p = parts[i].c_str();
            for (j = marker_num; j < marker_num + hap_len; j++) {
                psum->samples[sindex].nucs_1[j] = *p;
                p++;
            }
            i++;
            p = parts[i].c_str();
            for (j = marker_num; j < marker_num + hap_len; j++) {
                psum->samples[sindex].nucs_2[j] = *p;
                p++;
            }
            i++;
            sindex++;
        }

        marker_num += hap_len;

    } while (!gzeof(gz_fh));

    gzclose(gz_fh);

    return psum;
}

int
parse_population_map(string popmap_path, map<string, int> &pop_map, map<int, int> &pop_cnts)
{
    char   line[max_len];
    char   pop_id_str[id_len];
    vector<string> parts;
    uint   len;

    if (pmap_path.length() == 0)
        return 0;

    cerr << "Parsing population map.\n";

    ifstream fh(popmap_path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening population map '" << popmap_path << "'\n";
        return 0;
    }

    while (fh.good()) {
        fh.getline(line, max_len);

        len = strlen(line);
        if (len == 0) continue;

        //
        // Check that there is no carraige return in the buffer.
        //
        if (line[len - 1] == '\r') line[len - 1] = '\0';

        //
        // Ignore comments
        //
        if (line[0] == '#') continue;

        //
        // Parse the population map, we expect:
        // <file name> <tab> <population ID>
        //
        parse_tsv(line, parts);

        if (parts.size() != 2) {
            cerr << "Population map is not formated correctly: expecting two, tab separated columns, found " << parts.size() << ".\n";
            return 0;
        }

        strncpy(pop_id_str, parts[1].c_str(), id_len);
        for (int i = 0; i < id_len && pop_id_str[i] != '\0'; i++)
            if (!isdigit(pop_id_str[i])) {
                cerr << "Population map is not formated correctly: expecting numerical ID in second column, found '" << parts[1] << "'.\n";
                return 0;
            }

        //
        // Add the sample name to population number mapping.
        //
        pop_map[parts[0]] = atoi(parts[1].c_str());
        if (pop_cnts.count(atoi(parts[1].c_str())) == 0)
            pop_cnts[atoi(parts[1].c_str())] = 1;
        else
            pop_cnts[atoi(parts[1].c_str())]++;
    }

    fh.close();

    return 0;
}

int
build_file_list(vector<pair<int, string> > &files)
{
    vector<string> parts;
    string pattern;

    //
    // Read all the files from the Stacks directory.
    //
    uint   pos;
    string file;
    struct dirent *direntry;

    DIR *dir = opendir(in_path.c_str());

    if (dir == NULL) {
        cerr << "Unable to open directory '" << in_path << "' for reading.\n";
        exit(1);
    }

    switch(in_file_type) {
    case FileT::beagle:
        pattern = ".phased.gz";
        break;
    case FileT::fastphase:
    default:
        pattern = "_hapguess_switch.out";
        break;
    }

    while ((direntry = readdir(dir)) != NULL) {
        file = direntry->d_name;

        if (file == "." || file == "..")
            continue;

        pos = file.rfind(pattern);
        if (pos < file.length())
            files.push_back(make_pair(1, file.substr(0, pos)));
    }

    closedir(dir);

    if (files.size() == 0) {
        cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
        return 0;
    }

    return 1;
}

int parse_command_line(int argc, char* argv[]) {
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"haplotypes",  no_argument,       NULL, 'H'},
            {"skip-zeros",  no_argument,       NULL, 'Z'}, {"skip_zeros",  no_argument,       NULL, 'Z'},
            {"infile-type", required_argument, NULL, 't'}, {"infile_type", required_argument, NULL, 't'},
            {"num-threads", required_argument, NULL, 'p'}, {"num_threads", required_argument, NULL, 'p'},
            {"in-path",     required_argument, NULL, 'P'}, {"in_path",     required_argument, NULL, 'P'},
            {"cat-path",    required_argument, NULL, 'S'}, {"cat_path",    required_argument, NULL, 'S'},
            {"pop-map",     required_argument, NULL, 'M'}, {"pop_map",     required_argument, NULL, 'M'},
            {"batch-id",    required_argument, NULL, 'b'}, {"batch_id",    required_argument, NULL, 'b'},
            {"dprime-bin-size",   required_argument, NULL, 'B'}, {"dprime_bin_size",   required_argument, NULL, 'B'},
            {"minor-allele-freq", required_argument, NULL, 'a'}, {"minor_allele_freq", required_argument, NULL, 'a'},
            {"min-inform-pairs",  required_argument, NULL, 'm'}, {"min_inform_pairs",  required_argument, NULL, 'm'},
            {"dprime-threshold",  required_argument, NULL, 'T'}, {"dprime_threshold",  required_argument, NULL, 'T'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvZHAb:M:t:P:S:p:a:B:T:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'b':
            batch_id = is_integer(optarg);
            if (batch_id < 0) {
                cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
                help();
            }
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case 'a':
            minor_freq_lim = atof(optarg);
            break;
        case 'm':
            min_inform_pairs = atof(optarg);
            break;
        case 'P':
            in_path = optarg;
            break;
        case 'S':
            cat_path = optarg;
            break;
        case 't':
            if (strcasecmp(optarg, "phase") == 0)
                in_file_type = FileT::phase;
            else if (strcasecmp(optarg, "fastphase") == 0)
                in_file_type = FileT::fastphase;
            else if (strcasecmp(optarg, "beagle") == 0)
                in_file_type = FileT::beagle;
            else
                in_file_type = FileT::unknown;
            break;
        case 'M':
            pmap_path = optarg;
            break;
        case 'H':
            haplotypes = true;
            break;
        case 'Z':
            write_zeros = false;
            break;
        case 'B':
            bucket_dist = atoi(optarg);
            break;
        case 'T':
            dprime_threshold = true;
            dprime_threshold_level = atof(optarg);
            break;
        case 'v':
            version();
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;
        default:
            help();
            exit(1);
        }
    }

    if (in_path.length() == 0) {
        cerr << "You must specify a path to the directory containing Stacks output files.\n";
        help();
    }

    if (in_path.at(in_path.length() - 1) != '/')
        in_path += "/";

    if (minor_freq_lim > 0) {
        if (minor_freq_lim > 1)
            minor_freq_lim = minor_freq_lim / 100;

        if (minor_freq_lim > 0.5) {
            cerr << "Unable to parse the minor allele frequency\n";
            help();
        }
    }

    if (min_inform_pairs > 0) {
        if (min_inform_pairs > 1)
            min_inform_pairs = min_inform_pairs / 100;
    }

    return 0;
}

void version() {
    cerr << "phasedstacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "phasedstacks " << VERSION << "\n"
              << "phasedstacks -b id -S path -P path -t file_type [-p threads] [-M popmap] [-v] [-h]" << "\n"
              << "  b: Stacks batch ID.\n"
              << "  P: path to the phased output files.\n"
              << "  S: path to the Stacks output files.\n"
              << "  t: input file type. Supported types: fastphase, and beagle.\n"
              << "  p: number of processes to run in parallel sections of code.\n"
              << "  M: path to the population map, a tab separated file describing which individuals belong in which population.\n"
              << "  v: print program version." << "\n"
              << "  h: display this help messsage." << "\n"
              << "  --haplotypes: data were phased as RAD locus haplotypes.\n"
              << "  --dprime-bin-size: size of buckets for binning SNPs at a particular distance to calculate the mean D' value.\n"
              << "  --dprime-threshold <val>: if D' values fall above <val>, set the D' to 1, otherwise set D' to 0.\n\n"
              << "  Filtering options:\n"
              << "  --skip-zeros: do not include D' values of zero in the D' output.\n"
              << "  --minor-allele-freq: specify a minimum minor allele frequency required to process a nucleotide site (0 < a < 0.5).\n"
              << "  --min-inform-pairs: when building D' haplotype blocks, the minimum number of informative D' measures to combine two blocks (default 0.9).\n\n";


    exit(1);
}
