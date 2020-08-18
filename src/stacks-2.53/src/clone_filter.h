// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __CLONE_FILTER_H__
#define __CLONE_FILTER_H__

#include "constants.h"

#include <cstdlib>
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
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

#include "kmers.h"
#include "BustardI.h"   // Reading input files in Tab-separated Bustard format
#include "FastqI.h"     // Reading input files in FASTQ format
#include "FastaI.h"     // Reading input files in FASTA format
#include "gzFasta.h"    // Reading gzipped input files in FASTA format
#include "gzFastq.h"    // Reading gzipped input files in FASTQ format
#include "BamUnalignedI.h"
#include "utils.h"
#include "clean.h"
#include "file_io.h"
#include "write.h"

class Pair {
public:
    string p1_id;
    string p2_id;
    string p1_qual;
    string p2_qual;

    Pair(string p1_id, string p2_id, string p1_qual, string p2_qual) {
        this->p1_id   = p1_id;
        this->p2_id   = p2_id;
        this->p1_qual = p1_qual;
        this->p2_qual = p2_qual;
    }
    Pair(string p1_id, string p2_id) {
        this->p1_id   = p1_id;
        this->p2_id   = p2_id;
    }
};

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<char *, map<string, vector<Pair> >, hash_charptr, eqstr> CloneHash;
typedef sparse_hash_map<string, map<string, uint16_t> > OligoHash;
#else
typedef unordered_map<char *, map<string, vector<Pair> >, hash_charptr, eqstr> CloneHash;
typedef unordered_map<string, map<string, uint16_t> > OligoHash;
#endif

int  process_paired_reads(string, string, map<string, long> &, OligoHash &);
int  process_reads(string, map<string, long> &, OligoHash &);
int  process_paired_reads_by_sequence(string, string, map<string, long> &, CloneHash &, vector<char *> &);
int  write_clonereduced_sequence(string, string, CloneHash &, map<int, int> &, map<string, long> &);

int  free_hash(vector<char *> &);

void help( void );
void version( void );
int  parse_command_line(int, char**);

#endif // __CLONE_FILTER_H__
