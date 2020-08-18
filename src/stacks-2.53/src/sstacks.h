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

#ifndef __SSTACKS_H__
#define __SSTACKS_H__

#include "constants.h"

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <utility>

#include <string>

#include <iostream>
#include <fstream>

#include <vector>
#include <map>
#include <set>
#include <queue>
using std::queue;

#include <unordered_map>

#include "constants.h"
#include "kmers.h"
#include "stacks.h"
#include "locus.h"
#include "GappedAln.h"
#include "kmers.h"
#include "sql_utilities.h"
#include "utils.h"
#include "aln_utils.h"

typedef unordered_map<const char *, vector<pair<int, allele_type> >, hash_charptr, eqstr> HashMap;

void   help( void );
void   version( void );
int    parse_command_line(int, char**);
int    populate_hash(map<int, Locus *> &, HashMap &, vector<char *> &);
int    find_matches_by_sequence(map<int, Locus *> &, HashMap &, map<int, QLocus *> &);
int    find_matches_by_genomic_loc(map<int, Locus *> &, map<int, QLocus *> &);
int    verify_sequence_match(map<int, Locus *> &, QLocus *, set<int> &, map<string, vector<string> > &, unsigned long &, unsigned long &);
int    search_for_gaps(map<int, Locus *> &, map<int, QLocus *> &, KmerHashMap &, map<int, pair<allele_type, int> > &, double);
bool   verify_gapped_match(map<int, Locus *> &, QLocus *, set<int> &, map<allele_type, map<allele_type, AlignRes> > &, uint &, uint &, uint &, uint &, uint &);
int    verify_genomic_loc_match(Locus *, QLocus *, set<string> &, unsigned long &);
string generate_query_allele(Locus *, Locus *, const char *, allele_type);
bool   match_alleles(allele_type, allele_type);
int    generate_query_haplotypes(Locus *, QLocus *, set<string> &);
int    impute_haplotype(string, vector<pair<allele_type, string> > &, string &);
bool   compare_dist(pair<int, int>, pair<int, int>);
int    write_matches(string, map<int, QLocus *> &);

#endif // __SSTACKS_H__
