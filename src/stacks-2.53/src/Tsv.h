// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

//
// Code to parse the internal (and tempoary) data format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <chromosome> <base pair> <sequence> <phred quality score>
//
// One record per line.
//
#ifndef __TSV_H__
#define __TSV_H__

#include "input.h"

class Tsv: public Input {

 public:
    Tsv(const char *path) : Input(path) {};
    ~Tsv() {};
    Seq *next_seq();
    int  next_seq(Seq &);
};

Seq *Tsv::next_seq() {
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int Tsv::next_seq(Seq& s) {
    vector<string> parts;

    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
        return 0;
    }

    parse_tsv(this->line, parts);

    string id = parts[0] + "_" + parts[1];

    s = Seq(id.c_str(), parts[2].c_str(), parts[3].c_str(), parts[0].c_str(), atoi(parts[1].c_str()), strand_plus);

    return 1;
}

#endif // __TSV_H__
