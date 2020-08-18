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

#ifndef __CSTACKS_H__
#define __CSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif

#include <cerrno>
#include <zlib.h>   // Support for gzipped output files.

#include <getopt.h> // Process command-line options
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <utility>

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <queue>
using std::queue;
#include <algorithm>

#include "constants.h"
#include "stacks.h"
#include "kmers.h"
#include "locus.h"
#include "GappedAln.h"
#include "sql_utilities.h"
#include "aln_utils.h"
#include "utils.h"

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  initialize_new_catalog(pair<int, string> &, map<int, CLocus *> &);
int  initialize_existing_catalog(string, map<int, CLocus *> &);
int  update_catalog_index(map<int, CLocus *> &, map<string, int> &);
int  merge_catalog_loci(map<int, CLocus *> &, vector<int> &);
int  find_kmer_matches_by_sequence(map<int, CLocus *> &, map<int, QLocus *> &, int);
int  search_for_gaps(map<int, CLocus *> &, map<int, QLocus *> &, double, double);
int  find_matches_by_sequence(map<int, CLocus *> &, map<int, QLocus *> &);
int  find_matches_by_genomic_loc(map<string, int> &, map<int, QLocus *> &);
int  characterize_mismatch_snps(CLocus *, QLocus *);
int  merge_allele(Locus *, SNP *);
int  merge_matches(map<int, CLocus *> &, map<int, QLocus *> &, pair<int, string> &, int, uint &, uint &, uint &, uint &, uint &, uint &);
int  add_unique_tag(pair<int, string> &, map<int, CLocus *> &, QLocus *);
bool compare_dist(pair<int, int>, pair<int, int>);

int  write_catalog(map<int, CLocus *> &);
int  write_simple_output(CLocus *, ofstream &, ofstream &, ofstream &);
int  write_gzip_output(CLocus *, gzFile &, gzFile &, gzFile &);

bool compare_matches(Match *, Match *);

int  populate_kmer_hash(map<int, CLocus *> &, CatKmerHashMap &, vector<char *> &, int);

#endif // __CSTACKS_H__
