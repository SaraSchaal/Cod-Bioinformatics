// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2019, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __GENOTYPE_DICTIONARIES_H__
#define __GENOTYPE_DICTIONARIES_H__

#include <string>
#include <map>

#include "constants.h"

enum class CrossT  {unk, dh, cp, bc1, f2};
enum class FormatT {unk, rqtl, joinmap, onemap};

void initialize_dictionaries(map<string, map<string, string> > &global_dictionary);
void load_cp_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_joinmap_cp_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_onemap_cp_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_bc_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_f2_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_mm_bc_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_mm_f2_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_dh_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_mm_dh_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary);
void load_segregation_ratios(CrossT type, map<string, map<string, double> > &segregation_ratios);

extern int encoded_gtypes[4][4];

inline
int
encode_gtype(char a)
{
    switch (a) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    }
    return -1;
}

#endif // __GENOTYPE_DICTIONARIES_H__
