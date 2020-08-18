// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2017, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __GZFASTA_H__
#define __GZFASTA_H__

#include "input.h"

#ifdef HAVE_LIBZ

#include <cerrno>
#include <zlib.h>

class GzFasta: public Input {
    gzFile gz_fh;
    string buf;

 public:
    GzFasta(const char *path);
    GzFasta(string path) : GzFasta(path.c_str()) {}
    GzFasta() : Input() { this->gz_fh = NULL; }
    ~GzFasta() { gzclose(this->gz_fh); }
    void open(string path) { this->open(path.c_str()); }
    void open(const char *path);
    Seq *next_seq();
    int  next_seq(Seq &);
};

#else  // If HAVE_LIBZ is undefined and zlib library is not present.

#include "input.h"

class GzFasta: public Input {
 public:
    GzFasta(const char *path) : Input() { cerr << "Gzip support was not enabled when Stacks was compiled.\n"; };
    ~GzFasta() {};
    Seq *next_seq()      { return NULL; }
    int  next_seq(Seq &) { return 0; }
};

#endif // HAVE_LIBZ

#endif // __GZFASTA_H__
