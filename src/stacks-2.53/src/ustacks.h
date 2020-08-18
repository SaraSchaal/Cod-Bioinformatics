// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2018, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __USTACKS_H__
#define __USTACKS_H__

#include "constants.h"

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif

#include <getopt.h> // Process command-line options
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> // std::setprecision

#include <vector>
#include <map>
#include <unordered_map>
#include <queue>
using std::queue;
#include <set>
#include <unistd.h>

#include "config.h"
#include "kmers.h"
#include "utils.h"
#include "DNASeq.h"     // Class for storing two-bit compressed DNA sequences
#include "stacks.h"     // Major data structures for holding stacks
#include "mstack.h"
#include "mst.h"        // Minimum spanning tree implementation
#include "models.h"     // Contains maximum likelihood statistical models.
#include "FastaI.h"     // Reading input files in FASTA format
#include "FastqI.h"     // Reading input files in FASTQ format
#include "gzFasta.h"    // Reading gzipped input files in FASTA format
#include "gzFastq.h"    // Reading gzipped input files in FASTQ format
#include "aln_utils.h"
#include "GappedAln.h"

class HVal {
 public:
    vector<int> ids;

    int count() {
        return this->ids.size();
    }
    int add_id(int id) {
        this->ids.push_back(id);
        return 0;
    }
};

typedef unordered_map<DNANSeq, HVal> DNASeqHashMap;

void   help( void );
void   version( void );
int    parse_command_line(int, char**);
void   load_radtags(string, DNASeqHashMap &, size_t &);
int    load_seq_ids(vector<char *> &);
void   reduce_radtags(DNASeqHashMap &, map<int, Stack *> &, map<int, Rem *> &, size_t &, size_t &);
int    free_radtags_hash(DNASeqHashMap &, vector<DNANSeq *> &);
int    populate_merged_tags(map<int, Stack *> &, map<int, MergedStack *> &);
void   merge_stacks(map<int, MergedStack *> &, size_t &);
int    call_consensus(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &, bool);
int    call_alleles(MergedStack *, vector<DNANSeq *> &, vector<read_type> &);
int    update_consensus(MergedStack *, map<int, Stack *> &, map<int, Rem *> &);
size_t merge_remainders(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);
int    search_for_gapped_remainders(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);
size_t merge_gapped_remainders(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);
int    write_results(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);

//
// Match MergedStacks using a k-mer hashing algorithm
//
int  calc_kmer_distance(map<int, MergedStack *> &, int);
int  search_for_gaps(map<int, MergedStack *> &);
int  merge_gapped_alns(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &);
int  edit_gapped_seqs(map<int, Stack *> &, map<int, Rem *> &, MergedStack *, vector<pair<char, uint> > &);
int  edit_gaps(vector<pair<char, uint> > &, char *);
int  dist(MergedStack *, MergedStack *, vector<pair<char, uint> > &);
bool rank_alignments(Aln, Aln);
//
// Calculate depth of coverage statistics for stacks
//
void calc_coverage_distribution(map<int, MergedStack *> &, double &, double &, double &, double &);

//
// Dealing with lumberjack (huge) stacks
//
size_t remove_repetitive_stacks(map<int, MergedStack *> &);
int  deleverage(map<int, MergedStack *> &, set<int> &, int, vector<MergedStack *> &);

//
// Debugging
//
int  dump_unique_tags(map<int, Stack *> &);
int  dump_merged_tags(map<int, MergedStack *> &);
int  dump_stack_graph(string, map<int, Stack *> &, map<int, MergedStack *> &, vector<int> &, map<int, map<int, double> > &, map<int, set<int> > &);

//
// Utilities
//
MergedStack *merge_tags(MergedStack *, MergedStack *, int);
MergedStack *merge_tags(map<int, MergedStack *> &, set<int> &, int);
MergedStack *merge_tags(map<int, MergedStack *> &, int *, int, int);
long double factorial(int);

#endif // __USTACKS_H__
