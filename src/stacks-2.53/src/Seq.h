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

#ifndef SEQ_H
#define SEQ_H

#include "constants.h"

enum strand_type {strand_plus, strand_minus};
enum class AlnT {null, primary, secondary, supplementary};

class PhyLoc {
    char* chr_;
    static const char empty_str[1];
public:
    const char* chr() const {return chr_==NULL ? empty_str : chr_;}
    uint        bp;
    strand_type strand;

    PhyLoc() {
        chr_   = NULL;
        bp     = 0;
        strand = strand_plus;
    }
    PhyLoc(const PhyLoc& other)
        : chr_(NULL), bp(other.bp), strand(other.strand) {
        if (other.chr_ != NULL) {
            chr_ = new char[strlen(other.chr())+1];
            strcpy(chr_, other.chr_);
        }
    }
    PhyLoc& operator=(const PhyLoc& other) {PhyLoc cp (other); swap(*this, cp); return *this;}
    PhyLoc(const char *chr, uint bp) {
        this->chr_   = new char[strlen(chr)  + 1];
        this->bp     = bp;
        this->strand = strand_plus;
        strcpy(this->chr_,  chr);
    }
    PhyLoc(const char *chr, uint bp, strand_type strnd) {
        this->chr_   = new char[strlen(chr)  + 1];
        this->bp     = bp;
        this->strand = strnd;
        strcpy(this->chr_,  chr);
    }
    PhyLoc(const string& s); // Expects "CHROM:BP1(:[+-])?".

    ~PhyLoc() {
        if (chr_ != NULL)
            delete [] chr_;
    }

    void set(const char *ochr, uint obp, strand_type ostrand) {
        if (chr_ != NULL)
            delete[] chr_;
        chr_   = new char[strlen(ochr)+1];
        strcpy(chr_,  ochr);
        bp     = obp;
        strand = ostrand;
    }

    void clear() {if(chr_!=NULL) {delete[] chr_; chr_=NULL;} bp=0; strand=strand_plus;}
    bool empty() const {return chr_ == NULL || strlen(chr_) == 0;}
    friend void swap(PhyLoc& p, PhyLoc& q);
    bool operator==(const PhyLoc& other) const;
    bool operator<(const PhyLoc& other) const;
};

class Seq {
public:
    char *id;
    long capacity;
    char *seq;
    char *qual;
    string comment;

    //
    // Information for an aligned sequence.
    //
    AlnT   aln_type;
    double pct_clipped;
    int    map_qual;
    char  *loc_str;
    PhyLoc loc;

    Seq();
    Seq(const Seq& other);
    Seq(const char *, const char *);
    Seq(const char *, const char *, const char *);
    Seq(const char *, const char *, const char *, const char *, uint, strand_type);
    Seq(const char *, const char *, const char *, const char *, uint, strand_type, AlnT, double, int);
    ~Seq( void ) {
        if (id != NULL)
            delete[] id;
        if (seq != NULL)
            delete[] seq;
        if (qual != NULL)
            delete[] qual;
        if (loc_str != NULL)
            delete[] loc_str;
    }

    void reserve(size_t n, bool with_qual);

    friend void swap(Seq&, Seq&);
    Seq& operator=(Seq&& other) {swap(*this, other); return *this;}
    Seq& operator=(const Seq& other) = delete;

    // delete_seq(): Voids the `seq` and `qual` members.
    void delete_seq() {
        if (seq != NULL) {
            delete[] seq;
            seq = NULL;
        }
        if (qual != NULL) {
            delete[] qual;
            qual = NULL;
        }
    }
};

//
// Inline definitions
// ===========
//

inline
void Seq::reserve(size_t n, bool with_qual) {
    if (capacity < long(n)) {
        capacity = n;
        delete seq;
        seq = new char[capacity + 1];
        *seq = '\0';
        if (with_qual) {
            delete qual;
            qual = new char[capacity + 1];
            *qual = '\0';
        } else if (qual != NULL) {
            // Delete it to keep the capacity predictable.
            delete qual;
            qual = NULL;
        }
    }
}

inline
void swap(PhyLoc& p, PhyLoc& q) {
    char* chr = p.chr_;
    p.chr_ = q.chr_;
    q.chr_ = chr;

    const uint bp = p.bp;
    p.bp = q.bp;
    q.bp = bp;

    const strand_type strand = p.strand;
    p.strand = q.strand;
    q.strand = strand;
}

inline
bool PhyLoc::operator==(const PhyLoc& other) const {
    if (bp == other.bp
            && strand == other.strand
            && strcmp(chr(), other.chr()) == 0)
        return true;
    else
        return false;
}

inline
bool PhyLoc::operator<(const PhyLoc& other) const {
    const int chrcmp = strcmp(chr(), other.chr());
    if (chrcmp != 0)
        // Alphanumeric.
        return chrcmp < 0;
    else if (bp != other.bp)
        return bp < other.bp;
    else
        // Minus strand first.
        return strand == strand_minus && other.strand == strand_plus;
}

#endif
