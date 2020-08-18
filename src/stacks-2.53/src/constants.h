// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

//
// Pull in the configuration variables from the configure script
//
#ifndef HAVE_CONFIG_H
#error "Configuration incomplete. (HAVE_CONFIG_H is undefined.)"
#endif

#include "config.h"

#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cfloat>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <array>
#include <string>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <algorithm>
#include <functional>
#include <memory>
#include <regex>

// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Debugging-related macros.
#ifdef DEBUG
#define IF_NDEBUG_TRY
#define IF_NDEBUG_CATCH_ALL_EXCEPTIONS

#else
//#define NDEBUG // xxx Uncomment before releasing v2.1
#define IF_NDEBUG_TRY \
    try {
#define IF_NDEBUG_CATCH_ALL_EXCEPTIONS \
    } catch (const std::exception& e) { \
        std::cerr << "Aborted."; \
        if (typeid(e) != typeid(std::exception)) \
            std::cerr << " (" << e.what() << ")"; \
        std::cerr << "\n"; \
        return 13; \
    }
#endif //DEBUG

#define DOES_NOT_HAPPEN \
    do{cerr << "At " << __FILE__ << ":" << __LINE__ << " This should never happen.\n"; throw exception();} while(false) \
    // n.b. do{..}while(false) requests a trailing ';' ({..} doesn't).

using std::vector;
using std::array;
using std::string;
using std::set;
using std::map;
using std::unordered_set;
using std::unordered_map;
using std::ostream;
using std::istream;
using std::streambuf;
using std::cout;
using std::cerr;
using std::cin;
using std::flush;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::pair;
using std::make_pair;
using std::stoi;
using std::to_string;
using std::size_t;
using std::getline;
using std::exception;
using std::move;
using std::unique_ptr;
using std::ios;

typedef unsigned int uint;
typedef unsigned char uchar;

//
// Maximum line length for parsing input files.
//
const int max_len = 1024;

//
// Maximum length of idetifiers, such as sequence IDs and chromosome names.
//
const int id_len = 255;

//
// Size to use for internal buffer size for gzipped files being read with libz.
//
const int libz_buffer_size = 1048576;

//
// Number of digits to use when writing numbers in scientific notation.
//
const unsigned int fieldw = 5;

//
// Supported file types
//
enum class FileT {unknown,
    sql,     gzsql,
    fasta,   gzfasta,
    fastq,   gzfastq,
    bowtie,  sam, bam, tsv,
    bustard, phase, fastphase, beagle,
    vcf, gzvcf
};

int stacks_handle_exceptions(const exception& e);

string remove_suffix(FileT, const string&);

FileT guess_file_type(const string&);

void escape_char(char c, string& s);

#endif
