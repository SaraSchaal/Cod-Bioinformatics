// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __DNANSeq_H__
#define __DNANSeq_H__

#include <cstring>
#include <string>
#include <functional>

#ifdef HAVE_SPARSEHASH
#include <tr1/functional>
#endif

#include "constants.h"

/*
#define BITMASK(b)     (1 << ((b) % CHAR_BIT))
#define BITSLOT(b)     ((b) / CHAR_BIT)
#define BITSET(a, b)   ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b)  ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb)  ((nb + CHAR_BIT - 1) / CHAR_BIT)
*/

//
// DNANSeq
// Compressed DNA Sequence Class. Handles N's.
//
class DNANSeq {
    // Three-bit compression (2.67 bases/byte)
    //    A == 000
    //    C == 001
    //    G == 010
    //    T == 011
    //    N == 100
    const static uint bits_per_nuc = 3;

    // The number of bits over which the sequence is stored.
    uint bits;

    // The array of bits. Padding bits, if any, are 0.
    unsigned char *s;

public:
#ifdef HAVE_SPARSEHASH
    DNANSeq() : bits(3), s(new unsigned char[1]) { *s = (1<<2); } // "N"
#endif
    DNANSeq(uint len, const char* str);
    DNANSeq(const DNANSeq& other);
    ~DNANSeq() {delete[] s;}
    DNANSeq& operator=(const DNANSeq&);

    char operator[](uint pos) const;
    uint size() const {return bits / bits_per_nuc;}
    void seq(char* buf) const;
    string seq() const;

    void extend(uint n_before, uint n_after);

    bool operator==(const DNANSeq& other) const;
    bool operator<(const DNANSeq& other) const;
    friend struct std::hash<DNANSeq>;

private:
    uint nbytes() const {return (bits+7)/8;}

    static unsigned char testbit(const unsigned char* s_, uint index) { return s_[index/8] & (1 << index%8); }
    static void setbit(unsigned char* s_, uint index) { s_[index/8] |= (1 << index%8); }
};

inline
bool DNANSeq::operator== (const DNANSeq& other) const {
    if (bits != other.bits)
        return false;
    return memcmp(s, other.s, nbytes()) == 0;
}

inline
bool DNANSeq::operator<(const DNANSeq& other) const {
    if (bits != other.bits)
        return bits < other.bits;
    else
        return memcmp(s, other.s, nbytes()) < 0;
}

// Specialization for std::hash
// Based on GCC
namespace std {
template<>
struct hash<DNANSeq> {
    size_t operator()(const DNANSeq& seq) const {
        size_t __result = static_cast<size_t>(14695981039346656037ULL);
        for (uint i = 0; i < seq.nbytes(); i++) {
            __result ^= static_cast<size_t>(seq.s[i]);
            __result *= static_cast<size_t>(1099511628211ULL);
        }
        return __result;
    }
};
#ifdef HAVE_SPARSEHASH
namespace tr1 {
template<>
struct hash<DNANSeq> {
    size_t operator()(const DNANSeq& seq) const { return std::hash<DNANSeq>()(seq); }
};
}
#endif
}

#endif // __DNANSeq_H__
