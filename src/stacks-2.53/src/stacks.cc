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
// stacks.cc -- routines for the stack-holding containers
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include <cassert>
#include "constants.h"

#include "stacks.h"

Rem::Rem() {
    this->id         = 0;
    this->seq        = NULL;
    this->utilized   = false;
}

Rem::Rem(int id, uint seq_id, DNANSeq *seq) {
    this->id       = id;
    this->utilized = false;

    this->map.push_back(seq_id);

    this->seq = new DNANSeq(*seq);
}

int Rem::add_id(uint id) {
    this->map.push_back(id);

    return 0;
}

int Rem::add_seq(const DNANSeq *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(*seq);

    return 0;
}

int Rem::add_seq(const char *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(strlen(seq), seq);

    return 0;
}

int PStack::add_id(const char *id) {
    char *f = new char[strlen(id) + 1];
    strcpy(f, id);
    this->map.push_back(f);

    return 0;
}

int PStack::add_seq(const char *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(strlen(seq), seq);

    return 0;
}

int PStack::add_seq(const DNANSeq *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(*seq);

    return 0;
}

void PStack::extend(const PhyLoc& phyloc, int length) {
    if (this->seq == NULL)
        return;

    assert(strcmp(phyloc.chr(), loc.chr()) == 0
           && phyloc.strand == loc.strand);

    if (loc.strand == strand_plus) {
        assert(loc.bp >= phyloc.bp
               && loc.bp + seq->size() <= phyloc.bp + length);
        seq->extend(
                loc.bp - phyloc.bp,
                phyloc.bp + length - loc.bp - seq->size());
    } else {
        const int this_last = int(loc.bp) - int(seq->size()) + 1;
        const int target_last = int(phyloc.bp) - int(length) + 1;
        assert(loc.bp <= phyloc.bp
               && this_last >= target_last);
        seq->extend(
                phyloc.bp - loc.bp,
                this_last - target_last);
    }

    loc = phyloc;
}

int Stack::add_id(uint id) {
    this->map.push_back(id);

    return 0;
}

int Stack::add_seq(const char *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(strlen(seq), seq);

    return 0;
}

int Stack::add_seq(const DNANSeq *seq) {
    if (this->seq != NULL)
        delete this->seq;

    this->seq = new DNANSeq(*seq);

    return 0;
}