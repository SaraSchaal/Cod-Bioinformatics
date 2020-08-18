// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __INPUT_H__
#define __INPUT_H__

#include <cerrno>
#include <zlib.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "constants.h"
#include "stacks.h"

//
// The base class for all of our Input classes, such as Tsv, Fastq, Fasta, etc.
//
class Input {
 public:
    string   path;
    ifstream fh;
    char     line[max_len];

    Input();
    Input(const char *path);
    virtual ~Input();
    virtual Seq *next_seq() = 0;
    virtual int  next_seq(Seq &) = 0;
};

int   parse_tsv(const char *, vector<string> &);
int   parse_ssv(const char *, vector<string> &);
int   read_line(ifstream &, char **, int *);
int   read_gzip_line(gzFile &, char **, int *);
bool  is_comment(const char *);

inline
void strip_whitespace(std::string& s) {
    auto right = s.end();
    while (right != s.begin()) {
        --right;
        if (!std::isspace(*right)) {
            ++right;
            break;
        }
    }
    s.erase(right, s.end());

    auto left = s.begin();
    while (left != s.end() && std::isspace(*left))
        ++left;
    s.erase(s.begin(), left);
}

#endif // __INPUT_H__
