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
// sstacks -- search for occurances of stacks in a catalog of stacks.
//

#include <regex>

#include "constants.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "MetaPopInfo.h"
#include "gzFastq.h"
#include "BamI.h"

#include "sstacks.h"

using namespace std;

// Global variables to hold command-line options.
queue<string> samples;
string  catalog_path;
string  out_path;
FileT   in_file_type = FileT::sql;
int     num_threads  =  1;
int     samp_id      =  0;
bool    verify_haplotypes       = true;
bool    impute_haplotypes       = true;
bool    require_uniq_haplotypes = false;
bool    gapped_alignments       = true;
searcht search_type             = sequence;
bool    write_all_matches       = false;

double  min_match_len   = 0.80;
double  max_gaps        = 2.0;
int     gapped_kmer_len = 19;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    uint sample_cnt = samples.size();

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, Locus *>  catalog;
    bool compressed = false;
    int  res;

    cerr << "Searching for matches by sequence identity...\n";

    catalog_path += "catalog";
    res = load_loci(catalog_path, catalog, 0, false, compressed);

    if (res == 0) {
        cerr << "Error: Unable to parse catalog, '" << catalog_path << "'\n";
        throw exception();
    }

    KmerHashMap kmer_map;
    map<int, pair<allele_type, int> > allele_map;
    vector<char *> kmer_map_keys;

    //
    // Build a hash map out of the catalog, for exact matching.;
    //
    HashMap        catalog_exact_map;
    vector<char *> catalog_exact_map_keys;
    cerr << "Populating kmer dictionary for exact matches...";
    populate_hash(catalog, catalog_exact_map, catalog_exact_map_keys);
    cerr << "done.\n";

    if (gapped_alignments) {
        cerr << "Populating kmer dictionary for gapped alignments...";
        populate_kmer_hash(catalog, kmer_map, kmer_map_keys, allele_map, gapped_kmer_len);
        cerr << "done.\n";
    }

    string sample_path;
    int    i = 1;

    while (!samples.empty()) {
        map<int, QLocus *> sample;

        sample_path = samples.front();
        samples.pop();

        cerr << "\nProcessing sample '" << sample_path << "' [" << i << " of " << sample_cnt << "]\n";

        res = load_loci(sample_path, sample, 2, false, compressed);

        if (res == 0) {
            cerr << "Error: Unable to parse '" << sample_path << "'\n";
            throw exception();
        }

        in_file_type = compressed == true ? FileT::gzsql : FileT::sql;

        //
        // Assign the ID for this sample data.
        //
        samp_id = sample.begin()->second->sample_id;

        //dump_loci(catalog);
        //dump_loci(sample);

        cerr << "Searching for sequence matches...\n";
        find_matches_by_sequence(catalog, catalog_exact_map, sample);

        if (gapped_alignments) {
            cerr << "Searching for gapped alignments...\n";
            search_for_gaps(catalog, sample, kmer_map, allele_map, min_match_len);
        }

        write_matches(sample_path, sample);
        i++;

        //
        // Free memory associated with the sample.
        //
        for (map<int, QLocus *>::iterator j = sample.begin(); j != sample.end(); j++)
            delete j->second;
        sample.clear();
    }

    //
    // Free memory associated with the hash.
    //
    for (auto i = catalog_exact_map.begin(); i != catalog_exact_map.end(); i++)
        i->second.clear();
    catalog_exact_map.clear();
    for (uint i = 0; i < catalog_exact_map_keys.size(); i++)
        delete [] catalog_exact_map_keys[i];
    catalog_exact_map_keys.clear();

    if (gapped_alignments)
        free_kmer_hash(kmer_map, kmer_map_keys);

    //
    // Free memory associated with the catalog.
    //
    for (map<int, Locus *>::iterator j = catalog.begin(); j != catalog.end(); j++)
        delete j->second;

    cerr << "\nsstacks is done.\n";
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
find_matches_by_genomic_loc(map<int, Locus *> &sample_1, map<int, QLocus *> &sample_2)
{
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, Locus *>::iterator j;
    int  k;
    char id[id_len];

    //
    // Build a hash map out of the first sample (usually the catalog)
    //
    //
    // Create a map of the genomic locations of stacks in sample_1
    //
    cerr << "  Creating map of genomic locations...";

    map<string, set<int> > locations;
    for (j = sample_1.begin(); j != sample_1.end(); j++) {
        snprintf(id, id_len - 1, "%s|%d|%c",
                 j->second->loc.chr(),
                 j->second->loc.bp,
                 j->second->loc.strand == strand_plus ? '+' : '-');
        locations[id].insert(j->second->id);
    }

    cerr << "done.\n";

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample_2.begin(); i != sample_2.end(); i++)
        keys.push_back(i->first);

    //
    // Initialize some counters
    //
    unsigned long matches = 0;
    unsigned long nomatch = 0;
    unsigned long nosnps  = 0;
    unsigned long tot_hap = 0;
    unsigned long ver_hap = 0;

    #pragma omp parallel private(i, j, k, id)
    {
        unsigned long verified;
        #pragma omp for reduction(+:matches) reduction(+:tot_hap) reduction(+:ver_hap) reduction(+:nomatch) reduction(+:nosnps)
        for (k = 0; k < (int) keys.size(); k++) {

            i = sample_2.find(keys[k]);
            snprintf(id, id_len - 1, "%s|%d|%c",
                     i->second->loc.chr(),
                     i->second->loc.bp,
                     i->second->loc.strand == strand_plus ? '+' : '-');

            if (locations.count(id) > 0) {
                Locus *tag;
                set<int>::iterator loc_it;
                vector<pair<allele_type, string> >::iterator q;

                matches++;

                for (loc_it = locations[id].begin(); loc_it != locations[id].end(); loc_it++) {
                    tag = sample_1[*loc_it];

                    //
                    // Generate haplotypes for query tag relative to the catalog tag.
                    //
                    set<string> query_haplotypes;
                    generate_query_haplotypes(tag, i->second, query_haplotypes);
                    tot_hap += query_haplotypes.size() > 0 ? query_haplotypes.size() : 1;

                    if (verify_haplotypes) {
                        verified = verify_genomic_loc_match(tag, i->second, query_haplotypes, nosnps);
                        ver_hap += verified;
                        if (verified == 0) nomatch++;
                    } else {
                        i->second->add_match(tag->id, tag->strings.begin()->first);
                    }
                }
            }
        }
    }

    cerr << keys.size() << " sample loci matched against the catalog containing " << sample_1.size() << " loci.\n"
         << "  " << matches << " matching loci, " << nomatch << " contained no verified haplotypes.\n"
         << "  " << nosnps  << " loci contained SNPs unaccounted for in the catalog and were excluded.\n"
         << "  " << tot_hap << " total haplotypes examined from matching loci, " << ver_hap << " verified.\n";

    return 0;
}

int verify_genomic_loc_match(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes, unsigned long &nosnps) {
    vector<SNP *>::iterator i, j;

    //
    // We have found a match between the genomic location of s1 and s2. We now want
    // to verify that the haplotypes are consistent between the tags, i.e. they
    // have the same number and types of SNPs.
    //

    //
    // 1. First, if there are no SNPs present in either the query or catalog, just
    //    check that the strings match.
    //
    uint min_len = s1_tag->len > s2_tag->len ? s2_tag->len : s1_tag->len;

    if (s1_tag->snps.size() == 0 &&
        s2_tag->snps.size() == 0 &&
        strncmp(s1_tag->con, s2_tag->con, min_len) == 0) {
        s2_tag->add_match(s1_tag->id, "consensus");
        return 1;
    }

    //
    // 2. Second, we will check that the query locus (s2_tag) does not have any SNPs
    //    lacking in the catalog tag (s1_tag).
    //
    bool found;
    for (j = s2_tag->snps.begin(); j != s2_tag->snps.end(); j++) {
        found = false;
        //
        // SNP occurs in a column that is beyond the length of the catalog
        //
        if ((*j)->col > min_len - 1)
            continue;

        for (i = s1_tag->snps.begin(); i != s1_tag->snps.end(); i++) {
            if ((*i)->col == (*j)->col)
                found = true;
        }
        //
        // Query locus posses a SNP not present in the catalog.
        //
        if (found == false) {
            nosnps++;
            return 0;
        }
    }

    //
    // Finally, check that one of the constructed alleles matches the allele
    // passed in on the stack.
    //
    string cat_haplotype;
    vector<pair<allele_type, string> >::iterator c;
    set<string>::iterator a;

    uint matches = 0;
    for (a = query_haplotypes.begin(); a != query_haplotypes.end(); a++) {

        if (impute_haplotypes) {
            int res = impute_haplotype(*a, s1_tag->strings, cat_haplotype);

            if (res > 0) {
                //
                // If the matching haplotype was imputed, record the depths of the query alleles
                // under the new, imputed alleles.
                //
                if (s2_tag->alleles.count(cat_haplotype) == 0) {
                    if (s2_tag->alleles.count(*a) > 0)
                        s2_tag->alleles[cat_haplotype] = s2_tag->alleles[*a];
                    else
                        s2_tag->alleles[cat_haplotype] = s2_tag->depth;
                }
                //cerr << s2_tag->id << "; Adding cat haplotype: " << cat_haplotype << " based on depth of " << *a << ", " << s2_tag->alleles[cat_haplotype] << "\n";
                s2_tag->add_match(s1_tag->id, cat_haplotype);
                matches++;
            } else if (res < 0) {
                cerr << "  Failure imputing haplotype for catalog locus: " << s1_tag->id << " and query tag: " << s2_tag->id << "\n";
            }
        } else {
            for (c = s1_tag->strings.begin(); c != s1_tag->strings.end(); c++)
                if (*a == c->first) {
                    //cerr << "  Adding match between " << s1_tag->id << " and " << c->first << "\n";
                    s2_tag->add_match(s1_tag->id, c->first);
                    matches++;
                }
        }
    }

    return matches;
}

// int impute_haplotype(string query_haplotype,
//                   vector<pair<allele_type, string> > &cat_haplotypes,
//                   string &match) {

//     uint max_len = query_haplotype.length() > cat_haplotypes[0].first.length() ?
//      query_haplotype.length() :
//      cat_haplotypes[0].first.length();

//     //cerr << "Query len: " << query_haplotype.length() << "; Max length: " << max_len << "\n";

//     vector<string> cur, next;

//     for (uint i = 0; i < cat_haplotypes.size(); i++)
//      cur.push_back(cat_haplotypes[i].first);
//     match = "";

//     //
//     // Examine the haplotypes one SNP at a time. If we are able to uniquely
//     // determine the catalog haplotype that the query haplotype corresponds
//     // to, return it.
//     //
//     uint j = 0;
//     while (cur.size() > 1 && j < max_len) {

//      for (uint i = 0; i < cur.size(); i++) {
//          //cerr << "Comparing query[" << j << "]: '" << query_haplotype[j] << "' to catalog '" << cur[i][j] << "'\n";
//          if (query_haplotype[j] == cur[i][j]) {
//              //cerr << "  Keeping this haplotype.\n";
//              next.push_back(cur[i]);
//          }
//      }
//      cur = next;
//      next.clear();
//      j++;
//     }

//     //
//     // If there is only one left, make sure what we have of the haplotype does match
//     // and its not simply an erroneously called haplotype.
//     //
//     if (cur.size() == 1 &&
//      strncmp(cur[0].c_str(), query_haplotype.c_str(), max_len) == 0) {
//      match = cur[0];
//      return 1;
//     }

//     //
//     // If, after examining all the available SNPs in the query haplotype, there is
//     // still more than a single possible catalog haplotype, then we can't impute it.
//     //
//     return 0;
// }

int impute_haplotype(string query_haplotype,
                     vector<pair<allele_type, string> > &cat_haplotypes,
                     string &match) {

    if (cat_haplotypes.size() == 0) {
        cerr << "Warning: malformed catalog tag: missing haplotype information.\n";
        return -1;
    }

    //cerr << "Examining " << query_haplotype << "\n";
    uint max_len = query_haplotype.length() > cat_haplotypes[0].first.length() ?
        query_haplotype.length() :
        cat_haplotypes[0].first.length();

    //cerr << "Query len: " << query_haplotype.length() << "; Max length: " << max_len << "\n";

    vector<string> cur, next;
    uint match_cnt, no_n_cnt;

    for (uint i = 0; i < cat_haplotypes.size(); i++)
        cur.push_back(cat_haplotypes[i].first);
    match = "";

    //
    // Examine the haplotypes one SNP at a time. If we are able to uniquely
    // determine the catalog haplotype that the query haplotype corresponds
    // to, return it.
    //
    uint j = 0;
    while (cur.size() > 1 && j < max_len) {

        for (uint i = 0; i < cur.size(); i++) {
            //cerr << "Comparing query[" << j << "]: '" << query_haplotype << "' to catalog '" << cur[i] << "'\n";
            if (require_uniq_haplotypes && (query_haplotype[j] == cur[i][j] || query_haplotype[j] == 'N')) {
                //cerr << "  Keeping this haplotype.\n";
                next.push_back(cur[i]);
            } else if (query_haplotype[j] == cur[i][j]) {
                //cerr << "  Keeping this haplotype.\n";
                next.push_back(cur[i]);
            } //else {
                //cerr << "  Discarding this haplotype.\n";
            //}
        }
        cur = next;
        next.clear();
        j++;
    }

    //
    // If there is only one left, make sure what we have of the haplotype does match
    // and its not simply an erroneously called haplotype.
    //
    no_n_cnt  = 0;
    match_cnt = 0;
    if (cur.size() == 1) {
        if (require_uniq_haplotypes) {
            for (uint k = 0; k < max_len; k++)
                if (query_haplotype[k] != 'N') no_n_cnt++;
            for (uint k = 0; k < max_len; k++)
                if (cur[0][k] == query_haplotype[k]) match_cnt++;

            if (match_cnt == no_n_cnt) {
                //cerr << "Keeping " << query_haplotype << "\n";
                match = cur[0];
                return 1;
            }
        } else {
            if (strncmp(cur[0].c_str(), query_haplotype.c_str(), max_len) == 0) {
                match = cur[0];
                return 1;
            }
        }
    }

    //
    // If, after examining all the available SNPs in the query haplotype, there is
    // still more than a single possible catalog haplotype, then we can't impute it.
    //
    return 0;
}

int
generate_query_haplotypes(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes)
{
    //
    // Construct a set of haplotypes from the query locus relative to the catalog locus.
    // (The query locus already has a set of haplotypes, however, they don't necessarily
    //  account for all the SNPs in the catalog, so we will augment them with sequence
    //  from the consensus.)
    //
    if (s1_tag->snps.size() == 0 && s2_tag->snps.size() == 0)
        return 0;

    vector<pair<string, SNP *> >   merged_snps;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<pair<string, SNP *> >::iterator   k;
    vector<SNP *>::iterator i;

    for (i = s1_tag->snps.begin(); i != s1_tag->snps.end(); i++)
        columns[(*i)->col] = make_pair("catalog", *i);

    for (i = s2_tag->snps.begin(); i != s2_tag->snps.end(); i++) {
        //
        // Is this column already represented in the catalog?
        //
        if (columns.count((*i)->col))
            columns[(*i)->col] = make_pair("both", *i);
        else
            columns[(*i)->col] = make_pair("query", *i);
    }

    for (c = columns.begin(); c != columns.end(); c++)
        merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    map<string, int> converted_alleles;
    map<string, int>::iterator b;
    string old_allele, new_allele;
    int    pos;

    for (b = s2_tag->alleles.begin(); b != s2_tag->alleles.end(); b++) {
        old_allele = b->first;
        new_allele = "";
        pos        = 0;

        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            //
            // If the SNPs from the catalog haplotype beyond the length of the query, add Ns
            //
            if (k->first == "catalog") {
                new_allele += (k->second->col > s2_tag->len - 1) ? 'N' : s2_tag->con[k->second->col];
            } else {
                new_allele += old_allele[pos];
                pos++;
            }
        }
        query_haplotypes.insert(new_allele);
        converted_alleles[new_allele] = b->second;

        // cerr << "Adding haplotype: " << new_allele << " [" << b->first << "]\n";
    }

    if (s2_tag->alleles.size() == 0) {
        new_allele = "";
        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            new_allele += (k->second->col > s2_tag->len - 1) ? 'N' : s2_tag->con[k->second->col];
        }
        query_haplotypes.insert(new_allele);
        // cerr << "Adding haplotype 2: " << new_allele << "\n";
    } else {
        s2_tag->alleles.clear();
        for (b = converted_alleles.begin(); b != converted_alleles.end(); b++)
            s2_tag->alleles[b->first] = b->second;
    }

    return 0;
}

int
find_matches_by_sequence(map<int, Locus *> &sample_1, HashMap &sample_1_map, map<int, QLocus *> &sample_2)
{
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (auto i = sample_2.begin(); i != sample_2.end(); i++)
        keys.push_back(i->first);

    //
    // Initialize some counters
    //
    unsigned long matches = 0;
    unsigned long mmatch  = 0;
    unsigned long nosnps  = 0;
    unsigned long nomatch = 0;
    unsigned long no_haps = 0;
    unsigned long tot_hap = 0;
    unsigned long ver_hap = 0;

    #pragma omp parallel
    {
        #pragma omp for reduction(+:matches) reduction(+:tot_hap) reduction(+:ver_hap) reduction(+:nomatch) reduction(+:mmatch)
        for (uint k = 0; k < keys.size(); k++) {
            QLocus *query = sample_2[keys[k]];

            //
            // Iterate through the haplotypes for this tag in sample_2
            //
            HashMap::iterator hit;
            map<string, vector<string> > haplo_hits;
            set<int> loci_hit;

            for (auto q = query->strings.begin(); q != query->strings.end(); q++) {
                // cerr << "  Looking for haplotype: " << q->first << " with sequence " << q->second.substr(0, min_tag_len) << "\n";

                hit = sample_1_map.find(q->second.c_str());

                if (hit != sample_1_map.end()) {
                    tot_hap++;
                    // cerr << "    Found a match for " << hit->first << "\n";

                    for (auto c = hit->second.begin(); c != hit->second.end(); c++) {
                        //
                        // Record the catalog loci hit by the haplotypes of this query locus.
                        //
                        loci_hit.insert(c->first);

                        //
                        // Record the haplotypes hit between the query and catalog loci.
                        //
                        haplo_hits[q->first].push_back(c->second);

                        if (verify_haplotypes == false)
                            query->add_match(c->first, c->second);
                    }
                }
            }

            if (loci_hit.size() == 0)
                nomatch++;
            else if (loci_hit.size() > 0)
                matches++;

            if (verify_haplotypes && loci_hit.size() > 0) {
                uint verified = verify_sequence_match(sample_1, query, loci_hit, haplo_hits,
                                                      mmatch, nosnps);
                ver_hap += verified;
                if (verified == 0) {
                    no_haps++;
                    if (!loci_hit.empty())
                        assert(write_all_matches ^ query->matches.empty());
                }
            }
        }
    }

    cerr << keys.size() << " sample loci compared against the catalog containing " << sample_1.size() << " loci.\n"
         << "  " << matches << " matching loci, " << no_haps << " contained no verified haplotypes.\n"
         << "  " << mmatch  << " loci matched more than one catalog locus and were excluded.\n"
         << "  " << nosnps  << " loci contained SNPs unaccounted for in the catalog and were excluded.\n"
         << "  " << tot_hap << " total haplotypes examined from matching loci, " << ver_hap << " verified.\n";

    return 0;
}

int verify_sequence_match(map<int, Locus *> &sample_1, QLocus *query,
                          set<int> &loci_hit, map<string, vector<string> > &haplo_hits,
                          unsigned long &mmatch, unsigned long &nosnps) {
    //
    // 1. Check that this query locus matches just a single catalog locus.
    //
    if (loci_hit.size() > 1) {
        mmatch++;
        if (write_all_matches)
            for (int cloc_id : loci_hit)
                query->add_match(cloc_id, "multi");
        return 0;
    }

    Locus *cat = sample_1[*(loci_hit.begin())];

    //
    // 2. Make sure the query has no SNPs unaccounted for in the catalog.
    //
    bool found;

    for (auto i = query->snps.begin(); i != query->snps.end(); i++) {
        found = false;

        for (auto j = cat->snps.begin(); j != cat->snps.end(); j++) {
            if ((*i)->col == (*j)->col)
                found = true;
        }
        //
        // Query locus posses a SNP not present in the catalog.
        //
        if (found == false) {
            nosnps++;
            if (write_all_matches)
                query->add_match(cat->id, "extra_snp");
            return 0;
        }
    }

    //
    // 3. We want a one-to-one correspondance between a query haplotype and a
    //    catalog haplotype. This relationship fails when the catalog and query seqeunces
    //    are different lengths and the full length haplotype can not be determined.
    //
    char cigar[id_len];
    map<string, int> cat_hap, query_hap;

    for (auto it = haplo_hits.begin(); it != haplo_hits.end(); it++) {
        const string   &query_seq    = it->first;
        vector<string> &catalog_hits = it->second;

        query_hap[query_seq] = catalog_hits.size();
        for (uint j = 0; j < catalog_hits.size(); j++)
            cat_hap[catalog_hits[j]]++;
    }

    uint verified = 0;
    for (auto it = haplo_hits.begin(); it != haplo_hits.end(); it++) {
        const string   &query_seq    = it->first;
        vector<string> &catalog_hits = it->second;

        for (uint j = 0; j < catalog_hits.size(); j++) {
            if (cat_hap[catalog_hits[j]] == 1 &&
                query_hap[query_seq]     == 1) {
                verified++;

                // Create a CIGAR string for this simple match.
                snprintf(cigar, id_len - 1, "%uM", query->len);
                // Record the match.
                query->add_match(cat->id, catalog_hits[j], query_seq, 0, cigar);

                //
                // If the matching haplotype was imputed, record the depths of the query alleles
                // under the new, imputed alleles.
                //
                if (query->alleles.count(catalog_hits[j]) == 0) {
                    if (query->alleles.count(query_seq) > 0)
                        query->alleles[catalog_hits[j]] = query->alleles[query_seq];
                    else
                        query->alleles[catalog_hits[j]] = query->depth;
                }
            }
        }
    }

    if (verified == 0) {
        if (write_all_matches)
            query->add_match(cat->id, "none_verified");
        return 0;
    }

    return verified;
}

int
search_for_gaps(map<int, Locus *> &catalog, map<int, QLocus *> &sample,
                KmerHashMap &kmer_map, map<int, pair<allele_type, int> > &allele_map,
                double min_match_len)
{
    //
    // Search for loci that can be merged with a gapped alignment.
    //

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (auto it = sample.begin(); it != sample.end(); it++)
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(sample[keys[0]]->con);
    int num_kmers = con_len - gapped_kmer_len + 1;

    uint gapped_aln = 0;
    uint matches    = 0;
    uint mmatches   = 0;
    uint nomatches  = 0;
    uint ver_hap    = 0;
    uint tot_hap    = 0;
    uint bad_aln    = 0;
    uint no_haps    = 0;
    uint nosnps     = 0;

    #pragma omp parallel
    {
        QLocus                     *query;
        Locus                      *tag_2;
        KmerHashMap::iterator       h;
        AlignRes                    aln_res;
        vector<char *>              kmers;
        set<string>                 uniq_kmers;
        vector<int>                 hits;
        vector<pair<int, int>>      ordered_hits;
        uint                        hit_cnt, index, prev_id, allele_id, hits_size, stop, top_hit;
        pair<allele_type, int>      cat_hit;
        string                      query_allele, query_seq, cat_allele, cat_seq;
        map<string, vector<string>> haplo_hits;
        set<int>                    loci_hit;
        vector<pair<char, uint>>    cigar;

        GappedAln *aln = new GappedAln();

        initialize_kmers(gapped_kmer_len, num_kmers, kmers);

        #pragma omp for schedule(dynamic) reduction(+:matches) reduction(+:nomatches) reduction(+:mmatches) reduction(+:gapped_aln) reduction(+:ver_hap) \
                                          reduction(+:tot_hap) reduction(+:bad_aln) reduction(+:no_haps)
        for (uint i = 0; i < keys.size(); i++) {
            query = sample[keys[i]];

            //
            // If we already matched this locus to the catalog without using gapped alignments, skip it now.
            //
            if (query->matches.size() > 0)
                continue;

            gapped_aln++;

            map<allele_type, map<allele_type, AlignRes>> query_hits;

            loci_hit.clear();

            for (auto allele = query->strings.begin(); allele != query->strings.end(); allele++) {

                query_allele = allele->first;
                query_seq    = allele->second;
                tot_hap++;

                generate_kmers_lazily(allele->second.c_str(), gapped_kmer_len, num_kmers, kmers);

                //
                // We want to create a list of unique kmers to search with; otherwise, repetitive kmers will
                // generate, multiple, spurious hits in sequences with multiple copies of the same kmer.
                //
                uniq_kmers.clear();
                for (int j = 0; j < num_kmers; j++)
                    uniq_kmers.insert(kmers[j]);

                hits.clear();
                ordered_hits.clear();

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
                hits_size = hits.size();

                if (hits_size == 0)
                    continue;

                prev_id   = hits[0];
                index     = 0;

                do {
                    hit_cnt   = 0;
                    allele_id = prev_id;

                    while (index < hits_size && (uint) hits[index] == prev_id) {
                        hit_cnt++;
                        index++;
                    }

                    if (index < hits_size)
                        prev_id = hits[index];

                    ordered_hits.push_back(make_pair(allele_id, hit_cnt));

                } while (index < hits_size);

                if (ordered_hits.size() == 0)
                    continue;

                //
                // Process the hits from most kmer hits to least kmer hits.
                //
                sort(ordered_hits.begin(), ordered_hits.end(), compare_pair_intint);

                //
                // Only try to align the sequences with the most kmers in common.
                //
                top_hit = ordered_hits[0].second;
                stop    = 1;
                for (uint j = 1; j < ordered_hits.size(); j++)
                    if ((uint) ordered_hits[j].second < top_hit) {
                        stop = j;
                        break;
                    }

                for (uint j = 0; j < stop; j++) {
                    cat_hit = allele_map.at(ordered_hits[j].first);
                    hit_cnt = ordered_hits[j].second;

                    tag_2 = catalog[cat_hit.second];

                    cat_allele = cat_hit.first;
                    cat_seq    = "";
                    for (uint k = 0; k < tag_2->strings.size(); k++)
                        if (tag_2->strings[k].first == cat_hit.first) {
                            cat_seq = tag_2->strings[k].second;
                            break;
                        }

                    aln->init(tag_2->len, query->len);

                    // cerr << "Attempting to align: cat id " << tag_2->id << " with locus id " << query->id << "\n"
                    //      << "Cat allele: " << cat_allele   << "; seq: " << cat_seq << "\n"
                    //      << "Allele:     " << query_allele << "; seq: " << allele->second << "\n";

                    if (aln->align(cat_seq, query_seq)) {
                        aln->parse_cigar(cigar);

                        aln_res = aln->result();
                        //
                        // At this point in the analysis, all possible alleles we want to detect must already
                        // be present in the catalog. Therefore, we should reject any alignment that implies a
                        // change in the catalog sequence (with the presence of an deletion element) as
                        // spurious.
                        //
                        if (aln_res.cigar.find('D') != string::npos)
                            continue;

                        //
                        // If the alignment has too many gaps, skip it.
                        // If the alignment doesn't span enough of the two sequences, skip it.
                        //
                        if (aln_res.gap_cnt <= (max_gaps + 1) &&
                            aln_res.pct_id  >= min_match_len  &&
                            dist(cat_seq.c_str(), query_seq.c_str(), cigar) == 0) {
                            loci_hit.insert(tag_2->id);
                            query_hits[query_allele][cat_allele] = aln_res;
                        }
                    }
                }
            }

            if (verify_gapped_match(catalog, query, loci_hit, query_hits, mmatches, nosnps, no_haps, bad_aln, ver_hap))
                matches++;
            else {
                nomatches++;
                if (!loci_hit.empty())
                    assert(write_all_matches ^ query->matches.empty());
            }
        }

        //
        // Free the k-mers we generated for this query and the alignment class.
        //
        for (uint j = 0; j < kmers.size(); j++)
            delete [] kmers[j];
        kmers.clear();

        delete aln;
    }

    cerr << "Out of " << keys.size() << " query loci, " << gapped_aln << " gapped alignments attempted.\n"
         << "  "   << matches   << " loci matched one catalog locus; " << tot_hap << " total haplotypes examined, " << ver_hap << " verified.\n"
         << "  "   << nomatches << " loci matched no catalog locus;\n"
         << "    " << mmatches  << " loci matched more than one catalog locus and were excluded.\n"
         << "    " << nosnps    << " loci contained SNPs unaccounted for in the catalog and were excluded.\n"
         << "    " << no_haps   << " loci had no verified haplotypes.\n"
         << "    " << bad_aln   << " loci had inconsistent alignments to a catalog locus and were excluded.\n";

    return 0;
}

bool
verify_gapped_match(map<int, Locus *> &catalog, QLocus *query,
                    set<int> &loci_hit, map<allele_type, map<allele_type, AlignRes> > &query_hits,
                    uint &mmatch, uint &nosnps, uint &no_haps, uint &bad_aln, uint &ver_hits) {
    //
    // 1. Check that this query locus matches just a single catalog locus.
    //
    if (loci_hit.size() == 0) {
        return false;
    } else if (loci_hit.size() > 1) {
        mmatch++;
        if (write_all_matches)
            for (int cloc_id : loci_hit)
                query->add_match(cloc_id, "multi");
        return false;
    }

    int    cat_id = *(loci_hit.begin());
    Locus *cat    = catalog[cat_id];

    AlignRes aln_res;
    string   query_allele, cat_allele, converted_query_allele, qseq;
    uint     verified = 0;

    //
    // 2. We have aligned multiple query alleles against multiple catalog alleles for each locus.
    //    Different alleles may have slightly different CIGARs due to placement of gaps. Assign
    //    the best alignment for each query allele to a single catalog allele -- do not asign more
    //    than one query allele to more than one catalog allele.
    //
    set<string> cat_alleles_matched;
    Cigar cigar;
    for (auto query_it = query_hits.begin(); query_it != query_hits.end(); query_it++) {

        // Rank the alignments between this query allele and all catalog alleles.
        vector<pair<allele_type, AlignRes>> cat_hits;
        for (auto cat_it = query_it->second.begin(); cat_it != query_it->second.end(); cat_it++)
            cat_hits.push_back(*cat_it);
        sort(cat_hits.begin(), cat_hits.end(),
             [] (pair<allele_type, AlignRes> &a, pair<allele_type, AlignRes> &b)
             {
                 return compare_alignres(a.second, b.second);
             });

        for (uint i = 0; i < cat_hits.size(); i++)
            if (cat_alleles_matched.count(cat_hits[i].first) == 0) {
                query_it->second.clear();
                query_it->second.insert(cat_hits[i]);
                cat_alleles_matched.insert(cat_hits[i].first);
                break;
            }
    }

    //
    // 3. We need to apply the CIGAR to the query consensus sequence. Pick the match with the most depth
    //    so the consensus stays the consensus if there are different CIGARS between alleles.
    //
    string max_query_allele;
    int    max_query_allele_depth = 0;
    if (query->alleles.size() == 0) {
        max_query_allele = "consensus";
        max_query_allele_depth = query->depth;
    } else {
        for (auto query_it = query_hits.begin(); query_it != query_hits.end(); query_it++) {
            if (query->alleles[query_it->first] > max_query_allele_depth) {
                max_query_allele_depth = query->alleles[query_it->first];
                max_query_allele       = query_it->first;
            }
        }
    }
    parse_cigar(invert_cigar(query_hits[max_query_allele].begin()->second.cigar).c_str(), cigar);
    adjust_snps_for_gaps(cigar, query);
    qseq = apply_cigar_to_seq(query->con, cigar);
    query->add_consensus(qseq.c_str());

    //
    // 4. Make sure the query has no SNPs unaccounted for in the catalog.
    //
    bool found;
    for (auto i = query->snps.begin(); i != query->snps.end(); i++) {
        found = false;

        for (auto j = cat->snps.begin(); j != cat->snps.end(); j++) {
            if ((*i)->col == (*j)->col)
                found = true;
        }
        //
        // Query locus posses a SNP not present in the catalog.
        //
        if (found == false) {
            nosnps++;
            if (write_all_matches)
                query->add_match(cat->id, "extra_snp");
            return false;
        }
    }

    //
    // 5. Assign the allele hits after verifying there is a match between catalog and query allele.
    //
    for (auto query_it = query_hits.begin(); query_it != query_hits.end(); query_it++) {
        query_allele = query_it->first;
        cat_allele   = query_it->second.begin()->first;
        verified++;
        query->add_match(cat_id, cat_allele, query_allele, 0, invert_cigar(query_it->second.begin()->second.cigar));
    }

    if (verified > 0) {
        ver_hits += verified;
    } else {
        no_haps++;
        if (write_all_matches)
            query->add_match(cat->id, "none_verified");
        return false;
    }

    return true;
}

bool
match_alleles(allele_type catalog_allele, allele_type query_allele)
{
    const char *q    = catalog_allele.c_str();
    const char *p    = query_allele.c_str();
    const char *stop = p + query_allele.length();

    while (p < stop) {
        if (*p != 'N' && *p != *q)
            return false;
        p++;
        q++;
    }
    return true;
}

string
generate_query_allele(Locus *ctag, Locus *qtag, const char *qseq, allele_type allele)
{
    string new_allele = "";
    size_t qlen       = strlen(qseq);

    if (qtag->snps.size() == 0) {
        for (uint i = 0; i < ctag->snps.size(); i++)
            if (ctag->snps[i]->col > qlen - 1)
                new_allele += 'N';
            else if (qseq[ctag->snps[i]->col] == 'N')
                new_allele += ctag->con[ctag->snps[i]->col];
            else
                new_allele += qseq[ctag->snps[i]->col];
    } else {
        uint pos   = 0;
        uint index = 0;

        for (uint i = 0; i < ctag->snps.size(); i++) {
            if (index < qtag->snps.size() && qtag->snps[index]->col == ctag->snps[i]->col) {
                new_allele += allele[pos];
                index++;
                pos++;
            } else {
                if (ctag->snps[i]->col > qlen - 1)
                    new_allele += 'N';
                else if (qseq[ctag->snps[i]->col] == 'N')
                    new_allele += ctag->con[ctag->snps[i]->col];
                else
                    new_allele += qseq[ctag->snps[i]->col];
            }
        }
    }

    return new_allele;
}

int
populate_hash(map<int, Locus *> &sample, HashMap &hash_map, vector<char *> &hash_map_keys)
{
    Locus *tag;
    char  *key;

    //
    // Create a hash map out of the set of alleles for each Locus.
    //
    for (auto it = sample.begin(); it != sample.end(); it++) {
        tag = it->second;

        for (auto all_it = tag->strings.begin(); all_it != tag->strings.end(); all_it++) {
            key = new char[all_it->second.length() + 1];
            strncpy(key, all_it->second.c_str(), all_it->second.length());
            key[all_it->second.length()] = '\0';

            hash_map[key].push_back(make_pair(tag->id, all_it->first));
            hash_map_keys.push_back(key);
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
write_matches(string sample_path, map<int, QLocus *> &sample)
{
    map<int, QLocus *>::iterator i;

    //
    // Parse the input file names to create the output file
    //
    size_t pos_1    = sample_path.find_last_of("/");
    string out_file = out_path + sample_path.substr(pos_1 + 1)  + ".matches.tsv";

    if (in_file_type == FileT::gzsql)
        out_file += ".gz";

    //
    // Open the output files for writing.
    //
    gzFile   gz_matches=NULL;
    ofstream matches;
    if (in_file_type == FileT::gzsql) {
        gz_matches = gzopen(out_file.c_str(), "wb");
        if (!gz_matches) {
            cerr << "Error: Unable to open gzipped matches file '" << out_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_matches, libz_buffer_size);
        #endif
    } else {
        matches.open(out_file.c_str());
        check_open(matches, out_file);
    }

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
    log << "# sstacks version " << VERSION << "; generated on " << date << "\n";
    if (in_file_type == FileT::gzsql)
        gzputs_throwing(gz_matches, log.str().c_str());
    else
        matches << log.str();

    QLocus *qloc;
    string       type;
    uint         match_depth;
    stringstream sstr;

    cerr << "Outputing to file " << out_file << "\n";

    for (i = sample.begin(); i != sample.end(); i++) {
        qloc = i->second;

        for (uint j = 0; j < qloc->matches.size(); j++) {
            if (verify_haplotypes == false && search_type == genomic_loc)
                match_depth = qloc->depth;
            else
                match_depth =
                    qloc->alleles.count(qloc->matches[j]->cat_type) > 0 ?
                    qloc->alleles[qloc->matches[j]->cat_type] : qloc->depth;

            sstr << qloc->matches[j]->cat_id   << "\t"
                 << samp_id                    << "\t"
                 << qloc->id                   << "\t"
                 << qloc->matches[j]->cat_type << "\t"
                 << match_depth                << "\t"
                 << qloc->matches[j]->cigar    << "\n";
        }

        if (in_file_type == FileT::gzsql) gzputs_throwing(gz_matches, sstr.str().c_str()); else matches << sstr.str();
        sstr.str("");
    }

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log.str("");
    log << "# sstacks completed on " << date << "\n";
    if (in_file_type == FileT::gzsql)
        gzputs_throwing(gz_matches, log.str().c_str());
    else
        matches << log.str();

    if (in_file_type == FileT::gzsql)
        gzclose_throwing(gz_matches);
    else
        matches.close();

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    string in_dir;
    string popmap_path;

    while (1) {
        static struct option long_options[] = {
            {"help",              no_argument, NULL, 'h'},
            {"version",           no_argument, NULL, 'v'},
            {"aligned",           no_argument, NULL, 'g'},
            {"verify-hap",        no_argument, NULL, 'x'}, {"verify_hap",        no_argument, NULL, 'x'},
            {"uniq-haplotypes",   no_argument, NULL, 'u'}, {"uniq_haplotypes",   no_argument, NULL, 'u'},
            {"disable-gapped",    no_argument, NULL, 'G'}, {"disable_gapped",    no_argument, NULL, 'G'},
            {"threads",     required_argument, NULL, 'p'},
            {"catalog",     required_argument, NULL, 'c'},
            {"sample",      required_argument, NULL, 's'},
            {"out-path",    required_argument, NULL, 'o'}, {"out_path",    required_argument, NULL, 'o'},
            {"in-path",     required_argument, NULL, 'P'}, {"in_path",     required_argument, NULL, 'P'},
            {"popmap",      required_argument, NULL, 'M'},
            {"write-all-matches", no_argument, NULL, 2001},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "hgGxuvs:c:o:p:P:M:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case 's':
            samples.push(optarg);
            break;
        case 'g':
            search_type = genomic_loc;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'c':
            catalog_path = optarg;
            break;
        case 'x':
            verify_haplotypes = false;
            break;
        case 'u':
            require_uniq_haplotypes = true;
            break;
        case 'G':
            gapped_alignments = false;
            break;
        case 'P':
            in_dir = optarg;
            break;
        case 'M':
            popmap_path = optarg;
            break;
        case 'v':
            version();
            break;
        case 2001: //write-all-matches
            write_all_matches = true;
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

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    if (in_dir.empty() && catalog_path.empty()) {
        cerr << "Error: You must specify one of -P or -c.\n";
        help();
    } else if (
               ((!in_dir.empty() || !popmap_path.empty()) // One of -P, or -M
                && (!catalog_path.empty() || !samples.empty() || !out_path.empty())) // and one of -c, -s or -o
            ) {
        cerr << "Error: Please do not mix run modes (-P/-M or -c/-s/-o).\n";
        help();
    }

    if (!in_dir.empty()) {
        if (popmap_path.empty()) {
            cerr << "Error: Please specify some input samples (-M).\n";
            help();
        }


        if (in_dir.back() != '/')
            in_dir += "/";

        // Set `catalog_path`.
        catalog_path = in_dir;

        // Set `samples`.
        if (!popmap_path.empty()) {
            MetaPopInfo popmap;
            popmap.init_popmap(popmap_path);
            for (const Sample& s : popmap.samples())
                samples.push(in_dir + s.name);
        }

        // Set `out_path`.
        out_path = in_dir;

    } else if (!catalog_path.empty()) {
        if (samples.size() == 0) {
            cerr << "You must specify at least one sample file.\n";
            help();
        }

        if (catalog_path.back() != '/')
            catalog_path += "/";

        if (out_path.length() == 0)
            out_path = ".";

        if (out_path.back() != '/')
            out_path += "/";
    }

    return 0;
}

void version() {
    cerr << "sstacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "sstacks " << VERSION << "\n"
              << "sstacks -P dir -M popmap [-p n_threads]" << "\n"
              << "sstacks -c catalog_path -s sample_path [-s sample_path ...] -o path [-p n_threads]" << "\n"
              << "  -P,--in-path: path to the directory containing Stacks files.\n"
              << "  -M,--popmap: path to a population map file from which to take sample names.\n"
              << "  -s,--sample: filename prefix from which to load sample loci." << "\n"
              << "  -c,--catalog: path to the catalog." << "\n"
              << "  -p,--threads: enable parallel execution with n_threads threads.\n"
              << "  -o,--out-path: output path to write results." << "\n"
              << "  -x: don't verify haplotype of matching locus." << "\n"
              << "\n"
              << "Gapped assembly options:\n"
              << "  --disable-gapped: disable gapped alignments between stacks (default: enable gapped alignments).\n"
              ;

#ifdef DEBUG
    cerr << "\n"
            "Debug options:\n"
            "  --write-all-matches: Write blacklisted matches. (Compatibility with the web interface?)\n";
#endif

    exit(1);
}
