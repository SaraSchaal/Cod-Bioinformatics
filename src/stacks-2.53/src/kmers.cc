// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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
// kmers.cc -- routines to generate and hash K-mers
//
#include "kmers.h"

int determine_kmer_length(int read_len, int dist) {
    int kmer_len, span, min_matches;

    //
    // If distance allowed between sequences is 0, then k-mer length equals read length.
    //
    if (dist == 0)
        return read_len;

    //
    // Longer k-mer lengths will provide a smaller hash, with better key placement.
    // Increase the kmer_len until we start to miss hits at the given distance. Then
    // back the kmer_len off one unit to get the final value.
    //
    for (kmer_len = 7; kmer_len < read_len; kmer_len += 2) {
        span = (kmer_len * (dist + 1)) - 1;

        min_matches = read_len - span;

        if (min_matches <= 0) break;
    }

    if (kmer_len >= read_len) {
        cerr << "Unable to find a suitable k-mer length for matching.\n";
        exit(1);
    }

    kmer_len -= 2;

    return kmer_len;
}

int calc_min_kmer_matches(int kmer_len, int dist, int read_len, bool exit_err) {
    int span, min_matches;

    span = (kmer_len * (dist + 1)) - 1;

    min_matches = read_len - span;

    if (min_matches <= 0) {
        cerr <<
            "Warning: combination of k-mer length (" << kmer_len << ") and edit distance (" << dist << ") allows for " <<
            "sequences to be missed by the matching algorithm.\n";
    }

    if (min_matches <= 0 && exit_err)
        exit(1);
    else if (min_matches <= 0)
        min_matches = 1;

    return min_matches;
}

int
initialize_kmers(int kmer_len, int num_kmers, vector<char *> &kmers)
{
    char *kmer;

    for (int i = 0; i < num_kmers; i++) {
        kmer = new char[kmer_len + 1];
        kmers.push_back(kmer);
    }

    return 0;
}

int
generate_kmers_lazily(const char *seq, uint kmer_len, uint num_kmers, vector<char *> &kmers)
{
    char *kmer;
    const char *k = seq;

    if (num_kmers > kmers.size()) {
        int new_kmers = num_kmers - kmers.size();

        for (int i = 0; i < new_kmers; i++) {
            kmer = new char[kmer_len + 1];
            kmers.push_back(kmer);
        }
    }

    for (uint i = 0; i < num_kmers; i++) {
        kmer = kmers.at(i);
        strncpy(kmer, k, kmer_len);
        kmer[kmer_len] = '\0';
        k++;
    }

    return 0;
}

int
generate_kmers(const char *seq, int kmer_len, int num_kmers, vector<char *> &kmers)
{
    char *kmer;
    const char *k = seq;

    for (int i = 0; i < num_kmers; i++) {
        kmer = new char[kmer_len + 1];
        strncpy(kmer, k, kmer_len);
        kmer[kmer_len] = '\0';
        kmers.push_back(kmer);
        k++;
    }

    return 0;
}

int generate_permutations(map<int, char **> &pstrings, int width) {
    int   i, j, rem, div, num;
    char *p;
    //
    // Given a k-mer that allows wildcards -- 'N' characters, we need to generate all
    // possible k-mers. To do so, we will generate a range of numbers that we convert to
    // base 4, assuming that 0 = 'A', 1 = 'C', 2 = 'G', 3 = 'T'.
    //
    const int base = 4;
    int range      = (int) pow(4, width);

    //
    // Create an array of strings to hold the permuted nucleotides.
    //
    char **strings = new char * [range];
    for (i = 0; i < range; i++)
        strings[i] = new char[width + 1];

    for (i = 0; i < range; i++) {
        for (j = 0; j < width; j++)
            strings[i][j] = 'A';
        strings[i][width] = '\0';
    }

    for (i = 0; i < range; i++) {
        //
        // Convert this number to base 4
        //
        p   = strings[i]; p += width - 1;
        num = i;
        do {
            div = (int) floor(num / base);
            rem = num % base;

            switch(rem) {
            case 0:
                *p = 'A';
                break;
            case 1:
                *p = 'C';
                break;
            case 2:
                *p = 'G';
                break;
            case 3:
                *p = 'T';
                break;
            }
            num = div;
            p--;
        } while (div > 0);
    }

    pstrings[width] = strings;

    return 0;
}

int
populate_kmer_hash(map<int, MergedStack *> &merged, KmerHashMap &kmer_map, vector<char *> &kmer_map_keys, int kmer_len)
{
    map<int, MergedStack *>::iterator it;
    MergedStack    *tag;
    vector<char *>  kmers;
    bool            exists;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(merged.begin()->second->con) - kmer_len + 1;

    for (it = merged.begin(); it != merged.end(); it++) {
        tag = it->second;

        // Don't compute distances for masked tags
        if (tag->masked) continue;

        generate_kmers(tag->con, kmer_len, num_kmers, kmers);

        // Hash the kmers
        for (int j = 0; j < num_kmers; j++) {
            exists = kmer_map.count(kmers[j]) == 0 ? false : true;

            kmer_map[kmers[j]].push_back(tag->id);

            if (exists)
                delete [] kmers[j];
            else
                kmer_map_keys.push_back(kmers[j]);
        }
        kmers.clear();
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
populate_kmer_hash(map<int, Locus *> &catalog, CatKmerHashMap &kmer_map, vector<char *> &kmer_map_keys, int kmer_len)
{
    map<int, Locus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    vector<char *> kmers;
    Locus         *tag;
    char          *hash_key;
    bool           exists;
    int            num_kmers;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        num_kmers = strlen(tag->con) - kmer_len + 1;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            for (int j = 0; j < num_kmers; j++) {
                hash_key = kmers[j];
                exists   = kmer_map.count(hash_key) == 0 ? false : true;

                kmer_map[hash_key].push_back(make_pair(allele->first, tag->id));

                if (exists)
                    delete [] kmers[j];
                else
                    kmer_map_keys.push_back(hash_key);
            }
            kmers.clear();
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
populate_kmer_hash(map<int, Locus *> &catalog, KmerHashMap &kmer_map, vector<char *> &kmer_map_keys, map<int, pair<allele_type, int> > &allele_map, int kmer_len)
{
    map<int, Locus *>::iterator it;
    KmerHashMap::iterator   map_it;
    vector<pair<allele_type, string> >::iterator allele;
    map<int, pair<allele_type, int> >::iterator  allele_it;
    vector<char *> kmers;
    Locus         *tag;
    char          *hash_key;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers;
    int allele_index = 0;

    allele_it = allele_map.begin();

    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        num_kmers = strlen(tag->con) - kmer_len + 1;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            allele_it = allele_map.insert(allele_it, make_pair(allele_index, make_pair(allele->first, tag->id)));

            for (int j = 0; j < num_kmers; j++) {
                hash_key = kmers[j];

                map_it = kmer_map.find(hash_key);

                if (map_it != kmer_map.end()) {
                    map_it->second.push_back(allele_index);
                    delete [] kmers[j];
                } else {
                    kmer_map[hash_key].push_back(allele_index);
                    kmer_map_keys.push_back(hash_key);
                }
            }
            kmers.clear();

            allele_index++;
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
populate_kmer_hash(map<int, CLocus *> &catalog, KmerHashMap &kmer_map, vector<char *> &kmer_map_keys, map<int, pair<allele_type, int> > &allele_map, int kmer_len)
{
    map<int, CLocus *>::iterator it;
    KmerHashMap::iterator   map_it;
    vector<pair<allele_type, string> >::iterator allele;
    map<int, pair<allele_type, int> >::iterator  allele_it;
    vector<char *> kmers;
    Locus         *tag;
    char          *hash_key;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers;
    int allele_index = 0;

    allele_it = allele_map.begin();

    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        num_kmers = strlen(tag->con) - kmer_len + 1;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            allele_it = allele_map.insert(allele_it, make_pair(allele_index, make_pair(allele->first, tag->id)));

            for (int j = 0; j < num_kmers; j++) {
                hash_key = kmers[j];

                map_it = kmer_map.find(hash_key);

                if (map_it != kmer_map.end()) {
                    map_it->second.push_back(allele_index);
                    delete [] kmers[j];
                } else {
                    kmer_map[hash_key].push_back(allele_index);
                    kmer_map_keys.push_back(hash_key);
                }
            }
            kmers.clear();

            allele_index++;
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
free_kmer_hash(CatKmerHashMap &kmer_map, vector<char *> &kmer_map_keys)
{
    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        kmer_map[kmer_map_keys[i]].clear();
    }
    kmer_map.clear();

    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        delete [] kmer_map_keys[i];
    }
    kmer_map_keys.clear();

    return 0;
}

int
free_kmer_hash(KmerHashMap &kmer_map, vector<char *> &kmer_map_keys)
{
    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        kmer_map[kmer_map_keys[i]].clear();
    }
    kmer_map.clear();

    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        delete [] kmer_map_keys[i];
    }
    kmer_map_keys.clear();

    return 0;
}

int dist(const char *tag_1, Locus *tag_2, allele_type allele) {
    int   dist = 0;
    const char *p     = tag_1;
    const char *p_end = p + strlen(p);
    const char *q     = NULL;
    //
    // Identify which matching string has the proper allele
    //
    vector<pair<allele_type, string> >::iterator it;

    for (it = tag_2->strings.begin(); it != tag_2->strings.end(); it++)
        if (it->first == allele)
            q = it->second.c_str();
    if (q == NULL) return -1;

    const char *q_end = q + strlen(q);

    // Count the number of characters that are different
    // between the two sequences.
    while (p < p_end && q < q_end) {
        dist += (*p == *q) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

int
dist(const char *tag_1, const char *tag_2, vector<pair<char, uint> > &cigar)
{
    uint  size = cigar.size();
    char  op;
    uint  dist, len_1, len_2, pos_1, pos_2, stop;
    int   mismatches = 0;

    len_1 = strlen(tag_1);
    len_2 = strlen(tag_2);
    pos_1 = 0;
    pos_2 = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            //
            // A deletion has occured in tag_1 relative to tag_2.
            //
            pos_2 += dist;
            break;
        case 'I':
            //
            // An insertion has occured in tag_1 relative to tag_2.
            //
            pos_1 += dist;
            break;
        case 'M':
            stop = pos_1 + dist;
            while (pos_1 < stop && pos_1 < len_1 && pos_2 < len_2) {
                if (tag_1[pos_1] != 'N' && tag_2[pos_2] != 'N' && tag_1[pos_1] != tag_2[pos_2])
                    mismatches++;
                pos_1++;
                pos_2++;
            }
            break;
        default:
            break;
        }
    }

    return mismatches;
}

int
dist(Locus *tag_1, Locus *tag_2)
{
    int   dist  = 0;
    char *p     = tag_1->con;
    char *q     = tag_2->con;
    char *p_end = p + tag_1->len;
    char *q_end = q + tag_2->len;

    if (tag_1->len != tag_2->len) {
        if (tag_1->len < tag_2->len)
            dist += tag_2->len - tag_1->len;
        else if (tag_1->len > tag_2->len)
            dist += tag_1->len - tag_2->len;
    }

    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
        dist += (*p == *q) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

int
dist(MergedStack *tag_1, MergedStack *tag_2)
{
    int   dist  = 0;
    char *p     = tag_1->con;
    char *q     = tag_2->con;
    char *p_end = p + tag_1->len;
    char *q_end = q + tag_2->len;

    //
    // If the sequences are of different lengths, count the missing
    // nucleotides as mismatches.
    //
    if (tag_1->len != tag_2->len) {
        if (tag_1->len < tag_2->len)
            dist += tag_2->len - tag_1->len;
        else if (tag_1->len > tag_2->len)
            dist += tag_1->len - tag_2->len;
    }

    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
        dist += (*p == *q) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

int
dist(MergedStack *tag_1, char *seq)
{
    int   dist  = 0;
    char *p     = tag_1->con;
    char *q     = seq;
    uint  q_len = strlen(q);
    char *p_end = p + tag_1->len;
    char *q_end = q + q_len;

    //
    // If the sequences are of different lengths, count the missing
    // nucleotides as mismatches.
    //
    if (tag_1->len != q_len) {
        if (tag_1->len < q_len)
            dist += q_len - tag_1->len;
        else if (tag_1->len > q_len)
            dist += tag_1->len - q_len;
    }

    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
        dist += (*p == *q) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

bool compare_dist(pair<int, int> a, pair<int, int> b) {
    return (a.second < b.second);
}

int
check_frameshift(MergedStack *tag_1, MergedStack *tag_2, size_t mismatches)
{
    size_t cnt  = 0;
    size_t diff = 0;
    char const* p     = tag_1->con;
    char const* q     = tag_2->con;
    char const* p_end = p + tag_1->len - 1;
    char const* q_end = q + tag_2->len - 1;

    //
    // Set pointers to the common end of the sequences.
    //
    if (tag_1->len != tag_2->len) {
        if (tag_1->len < tag_2->len) {
            diff = tag_2->len - tag_1->len;
            p_end -= diff;
        } else if (tag_1->len > tag_2->len) {
            diff   = tag_1->len - tag_2->len;
            q_end -= diff;
        }
    }

    //
    // Count the number of characters that are different
    // at the 3' end of the sequence to test for possible frameshifts.
    //
    size_t i = diff;
    while (p_end >= p && q_end >= q && i < mismatches) {
        cnt += (*p_end != *q_end) ? 1 : 0;
        p_end--;
        q_end--;
        i++;
    }

    return cnt;
}

int
check_frameshift(const char *tag_1, Locus *tag_2, allele_type allele, size_t mismatches)
{
    size_t cnt = 0;
    const char *p     = tag_1;
    const char *q     = NULL;
    //
    // Identify which matching string has the proper allele
    //
    vector<pair<allele_type, string> >::iterator it;

    for (it = tag_2->strings.begin(); it != tag_2->strings.end(); it++)
        if (it->first == allele)
            q = it->second.c_str();
    if (q == NULL) return -1;

    const char *p_end = p + strlen(p) - 1;
    const char *q_end = q + strlen(q) - 1;

    //
    // Count the number of characters that are different
    // at the 3' end of the sequence to test for possible frameshifts.
    //
    size_t i = 0;
    while (p_end >= p && q_end >= q && i < mismatches) {
        cnt += (*p_end != *q_end) ? 1 : 0;
        p_end--;
        q_end--;
        i++;
    }

    return cnt;
}

int
check_frameshift(MergedStack *tag_1, const char *seq, size_t mismatches)
{
    int   cnt     = 0;
    const char *p = tag_1->con;
    const char *q = seq;
    uint  q_len   = strlen(q);
    const char *p_end = p + tag_1->len - 1;
    const char *q_end = q + q_len - 1;

    //
    // If the sequences are of different lengths, count the missing
    // nucleotides as mismatches.
    //
    if (tag_1->len != q_len) {
        if (tag_1->len < q_len)
            q_end -= q_len - tag_1->len;
        else if (tag_1->len > q_len)
            p_end -= tag_1->len - q_len;
    }

    //
    // Count the number of characters that are different
    // at the 3' end of the sequence to test for possible frameshifts.
    //
    size_t i = 0;
    while (p_end >= p && q_end >= q && i < mismatches) {
        cnt += (*p_end != *q_end) ? 1 : 0;
        p_end--;
        q_end--;
        i++;
    }

    return cnt;
}

int dump_kmer_map(KmerHashMap &kmer_map) {
    KmerHashMap::iterator kit;
    vector<int>::iterator vit;

    cerr << kmer_map.size() << " keys in the map.\n";

    int i = 1;
    for (kit = kmer_map.begin(); kit != kmer_map.end(); kit++) {
        cerr << "Key #" << i << " " << kit->first << ": ";
        for (vit = (kit->second).begin(); vit != (kit->second).end(); vit++)
            cerr << " " << *vit;
        cerr << "\n";
        i++;

        if (i > 1000) break;
    }

    return 0;
}
