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

#ifndef __BOWTIEI_H__
#define __BOWTIEI_H__

//
// Code to parse Bowtie's alignment format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> <sequence> <phred quality score> <flag> <mismatches>
//
// One record per line.
//

#include "input.h"

class Bowtie: public Input {

 public:
    Bowtie(const char *path) : Input(path) {};
    ~Bowtie() {};
    Seq *next_seq();
    int  next_seq(Seq &);
};

Seq *Bowtie::next_seq() {
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int Bowtie::next_seq(Seq& s) {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
        return 0;
    }

    parse_tsv(this->line, parts);

    strand_type strand = parts[1] == "+" ? strand_plus : strand_minus;

    //
    // If the read was aligned on the reverse strand (and is therefore reverse complemented)
    // alter the start point of the alignment to reflect the right-side of the read, at the
    // end of the RAD cut site.
    //
    int bp = strand == strand_plus ? atoi(parts[3].c_str()) : atoi(parts[3].c_str()) + parts[4].length();

    s = Seq(parts[0].c_str(), parts[4].c_str(), parts[5].c_str(),
                     parts[2].c_str(), bp, strand);

    return 1;
}

#endif // __BOWTIEI_H__
