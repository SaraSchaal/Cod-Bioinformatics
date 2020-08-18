// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2016-2018, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __ALN_UTILS_H__
#define __ALN_UTILS_H__

#include <tuple>

#include "constants.h"

typedef vector<pair<char, uint>> Cigar;

extern const bool is_cigar_char[256];

ostream& operator<< (ostream&, const Cigar&);

string invert_cigar(string);
int    invert_cigar(Cigar &cigar);
int    convert_local_cigar_to_global(Cigar &cigar);
int    parse_cigar(const char*, Cigar&, bool check_correctness = false);
string apply_cigar_to_seq(const char*, const Cigar&);
string remove_cigar_from_seq(const char*, const Cigar&);
string apply_cigar_to_model_seq(const char*, const Cigar&);
int    apply_cigar_to_seq(char*, uint, const char*, Cigar&);
int    apply_cigar_to_model_seq(char*, uint, const char*, Cigar&);

std::tuple<uint,uint,uint> cigar_lengths(const Cigar&);
inline uint cigar_length_padded(const Cigar& c) {return std::get<0>(cigar_lengths(c));}
inline uint cigar_length_ref(const Cigar& c) {return std::get<1>(cigar_lengths(c));}
inline uint cigar_length_query(const Cigar& c) {return std::get<2>(cigar_lengths(c));}
void cigar_extend_right(Cigar&, size_t);
void cigar_extend_left(Cigar&, size_t);
void cigar_trim_left(Cigar&, size_t);

void cigar_simplify_to_MDI(Cigar&); // Makes all operations to be one of 'M', 'D' or 'I'.
void cigar_canonicalize_MDI_order(Cigar&); // Makes the cigar match [I][D][M[D][I]]+.
inline bool cigar_is_MDI(const Cigar& cig) {for(auto& op : cig) if(op.first!='M' && op.first!='D' && op.first!='I') return false; return true;}

//deprecated.
inline void simplify_cigar_to_MDI(Cigar& c) {cigar_simplify_to_MDI(c);}

//
// Inline definitions.
// ==========
//

inline
void cigar_extend_right(Cigar& cig, size_t len) {
    if (len == 0)
        return;
    if (cig.back().first == 'D')
        cig.back().second += len;
    else
        cig.push_back({'D', len});
}

inline
void cigar_extend_left(Cigar& cig, size_t len) {
    if (len == 0)
        return;
    if (cig.front().first == 'D')
        cig.front().second += len;
    else
        cig.insert(cig.begin(), {'D', len});
}

inline
void cigar_trim_left(Cigar& cig, size_t len) {
    assert(!cig.empty());
    assert(cigar_is_MDI(cig));

    if (len == 0) {
        return;
    } else if (cig.front().first == 'D' && cig.front().second >= len) {
        cig.front().second -= len;
        if (cig.front().second == 0)
            cig.erase(cig.begin());
        return;
    }
    // Trimming of the CIGAR isn't trivial.
    // We must remove as many M/D operations as necessary, and merge the I operations.
    // In practice, we advance an iterator until we've trimmed enough bases, while
    // recording seen I operations (if any), then we erase all operations until that
    // iterator position.
    Cigar::iterator op = cig.begin();
    assert(op != cig.end());
    size_t n_i = 0;
    size_t to_trim = len;
    if (op->first == 'I') {
        n_i += op->second;
        ++op;
        assert(op != cig.end());
    }
    while (to_trim > 0) {
        assert(op != cig.end());
        assert(op->first == 'M' || op->first == 'D');
        if (op->second > to_trim) {
            if (op->first == 'M')
                n_i += to_trim;
            op->second -= to_trim;
            break;
        } else {
            if (op->first == 'M')
                n_i += op->second;
            to_trim -= op->second;
        }

        ++op;
        if (op == cig.end()) {
            assert(len == cigar_length_ref(cig));
        } else if (op->first == 'I') {
            n_i += op->second;
            ++op;
        }
    }
    assert(n_i > 0); // `n_i == 0` would only happen for the trivial cases treated above.

    if (op != cig.begin()) {
        cig.front() = {'I', n_i};
        cig.erase(cig.begin() + 1, op);
    } else {
        cig.insert(cig.begin(), {'I', n_i});
    }
}

#endif  // __ALN_UTILS_H__
