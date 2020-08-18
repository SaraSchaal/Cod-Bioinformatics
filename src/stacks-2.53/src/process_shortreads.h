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

#ifndef __PROCESS_SHORTREADS_H__
#define __PROCESS_SHORTREADS_H__

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

#include "constants.h"
#include "clean.h"
#include "file_io.h"
#include "utils.h"
#include "log_utils.h"
#include "write.h"
#include "BustardI.h"       // Reading input files in Tab-separated Bustard format
#include "FastqI.h"         // Reading input files in FASTQ format
#include "gzFastq.h"        // Reading gzipped input files in FASTQ format
#include "BamUnalignedI.h"  // Reading data from unaligned BAM files

void help( void );
void version( void );
int  parse_command_line(int, char **);
template<typename fhType>
int  process_reads(string,
                   set<string> &, set<string> &,
                   map<BarcodePair, fhType *> &,
                   map<string, long> &, map<BarcodePair, map<string, long> > &);
template<typename fhType>
int  process_paired_reads(string, string,
                          set<string> &, set<string> &,
                          map<BarcodePair, fhType *> &,
                          map<BarcodePair, fhType *> &,
                          map<BarcodePair, fhType *> &,
                          map<BarcodePair, fhType *> &,
                          map<string, long> &, map<BarcodePair, map<string, long> > &);
int  process_singlet(RawRead *,
                     bool,
                     map<string, long> &, map<string, long> &);

int  dist(const char *, char *);
int  print_results(int, char **, vector<BarcodePair> &, map<string, map<string, long> > &, map<BarcodePair, map<string, long> > &);
int  compare_barcodes(pair<BarcodePair, int>, pair<BarcodePair, int>);

#endif // __PROCESS_SHORTREADS_H__
