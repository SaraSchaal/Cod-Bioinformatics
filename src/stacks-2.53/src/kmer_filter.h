// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2014, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __KMER_FILTER_H__
#define __KMER_FILTER_H__

#include "constants.h"

#include <cstdlib>
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <unordered_map>

#ifdef HAVE_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

#include "clean.h"
#include "utils.h"
#include "kmers.h"
#include "write.h"
#include "BustardI.h"   // Reading input files in Tab-separated Bustard format
#include "FastaI.h"     // Reading input files in FASTA format
#include "FastqI.h"     // Reading input files in FASTQ format
#include "gzFasta.h"    // Reading gzipped input files in FASTA format
#include "gzFastq.h"    // Reading gzipped input files in FASTQ format

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<char *, long, hash_charptr, eqstr> SeqKmerHash;
#else
typedef unordered_map<char *, long, hash_charptr, eqstr> SeqKmerHash;
#endif

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  build_file_list(vector<string> &, vector<pair<string, string> > &);
int  process_reads(string, string, SeqKmerHash &, map<string, long> &);
int  process_paired_reads(string, string, string, string, SeqKmerHash &, map<string, long> &);
int  print_results(map<string, map<string, long> > &);

//
// Functions to normalize read depth
//
int  normalize_reads(string, string, SeqKmerHash &, vector<char *> &, map<string, long> &);
int  normalize_paired_reads(string, string, string, string, SeqKmerHash &, vector<char *> &, map<string, long> &);
bool normalize_kmer_lookup(SeqKmerHash &, char *, char *, int, vector<char *> &);

//
// Functions for finding and removing reads with rare kmers
//
int  populate_kmers(vector<pair<string, string> > &, vector<pair<string, string> > &, SeqKmerHash &, vector<char *> &);
int  process_file_kmers(string, SeqKmerHash &, vector<char *> &);
int  generate_kmer_dist(SeqKmerHash &);
int  calc_kmer_median(SeqKmerHash &, double &, double &);
int  kmer_map_cmp(pair<char *, long>, pair<char *, long>);
int  kmer_lookup(SeqKmerHash &, char *, char *, int, int &, int &);
int  free_kmer_hash(SeqKmerHash &, vector<char *> &);

int  read_kmer_freq(string, SeqKmerHash &, vector<char *> &);
int  write_kmer_freq(string, SeqKmerHash &);

#endif // __KMER_FILTER_H__
