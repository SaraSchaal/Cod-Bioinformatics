// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
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
// Code to parse Illumina's Bustard file format. It takes the tab-separated form:
//
//  <machine> <run> <lane> <tile> <x> <y> <index> <read> <seq> <phred> <filter>
//
// One record per line.
//
#ifndef __BUSTARDI_H__
#define __BUSTARDI_H__

#include "input.h"

const uint num_bustd_fields = 11;

class Bustard: public Input {

 public:
    Bustard(const char *path) : Input(path) {};
    Bustard(string path) : Input(path.c_str()) {};
    ~Bustard() {};
    Seq *next_seq();
    int  next_seq(Seq &);
};

Seq *Bustard::next_seq() {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
        return NULL;
    }

    parse_tsv(this->line, parts);

    if (parts.size() != num_bustd_fields) {
        cerr << "Error parsing '" << this->path.c_str() << "' found " << parts.size() << " fields, but expecting " << num_bustd_fields << "). "
             << "Perhaps you should specify the input file type (-i)?\n";
        return NULL;
    }

    Seq *s = new Seq;
    s->seq = new char[parts[8].length() + 1];
    strcpy(s->seq,  parts[8].c_str());
    s->qual = new char[parts[9].length() + 1];
    strcpy(s->qual, parts[9].c_str());

    s->id = new char[id_len];
    sprintf(s->id, "@%s:%s:%s:%s:%s#%s/%s",
            parts[0].c_str(),
            parts[2].c_str(),
            parts[3].c_str(),
            parts[4].c_str(),
            parts[5].c_str(),
            parts[6].c_str(),
            parts[7].c_str());

    return s;
}

int Bustard::next_seq(Seq &s) {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
        return 0;
    }

    parse_tsv(this->line, parts);

    strcpy(s.seq,  parts[2].c_str());
    strcpy(s.qual, parts[3].c_str());

    sprintf(s.id, "@%s:%s:%s:%s:%s#%s/%s",
            parts[0].c_str(),
            parts[1].c_str(),
            parts[2].c_str(),
            parts[3].c_str(),
            parts[4].c_str(),
            parts[5].c_str(),
            parts[6].c_str());

    return 1;
}

#endif // __BUSTARDI_H__
