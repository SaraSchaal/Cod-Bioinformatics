// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __FILE_IO_H__
#define __FILE_IO_H__

#include <cstdlib>
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <cerrno>

#include "constants.h"
#include "clean.h"
#include "input.h"

//
// Command line options defined in process_radtags and process_shortreads.
//
extern FileT     in_file_type;
extern FileT     out_file_type;
extern barcodet  barcode_type;
extern bool      paired;
extern bool      interleaved;
extern bool      merge;
extern string    out_path;
extern string    in_file;
extern string    in_file_p1;
extern string    in_file_p2;
extern string    in_path_1;
extern string    in_path_2;

//
// Defined externally in process_radtags and process_shortreads.
//
void help( void );

int    build_file_list(vector<pair<string, string> > &);
int    load_barcodes(string, vector<BarcodePair> &,
                     set<string> &, set<string> &,
                     uint &, uint &, uint &, uint &);
int    open_files(vector<pair<string, string> > &,
                  vector<BarcodePair> &,
                  map<BarcodePair, ofstream *> &,
                  map<BarcodePair, ofstream *> &,
                  map<BarcodePair, ofstream *> &,
                  map<BarcodePair, ofstream *> &,
                  map<string, map<string, long> > &);
int    open_files(vector<pair<string, string> > &,
                  vector<BarcodePair> &,
                  map<BarcodePair, gzFile *> &,
                  map<BarcodePair, gzFile *> &,
                  map<BarcodePair, gzFile *> &,
                  map<BarcodePair, gzFile *> &,
                  map<string, map<string, long> > &);
int    close_file_handles(map<BarcodePair, ofstream *> &);
int    close_file_handles(map<BarcodePair, gzFile *> &);

#endif // __FILE_IO_H__
