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

#ifndef __STACKS_H__
#define __STACKS_H__

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include<iostream>
#include<sstream>

#include "constants.h"
#include "nucleotides.h"
#include "utils.h"
#include "Seq.h"
#include "DNASeq.h"
#include "DNANSeq.h"
#include "DNASeq4.h"
#include "MetaPopInfo.h"

typedef unsigned int uint;
typedef string allele_type;

enum snp_type    {snp_type_het, snp_type_hom, snp_type_unk, snp_type_discarded};
enum read_type   {primary, secondary};
enum searcht     {sequence, genomic_loc};
enum corr_type   {p_value, bonferroni_win, bonferroni_gen, no_correction};

class SNP {
 public:
    snp_type type;   // Heterozygous or homozygous
    uint     col;
    float    lratio;
    char     rank_1;
    char     rank_2;
    char     rank_3;
    char     rank_4;

    SNP() {
        col    = 0;
        lratio = 0.0;
        rank_1 = 0;
        rank_2 = 0;
        rank_3 = 0;
        rank_4 = 0;
    }
    SNP(const SNP &other) {
        col    = other.col;
        type   = other.type;
        lratio = other.lratio;
        rank_1 = other.rank_1;
        rank_2 = other.rank_2;
        rank_3 = other.rank_3;
        rank_4 = other.rank_4;
    }
};

inline
bool compare_pair_snp(const pair<string, SNP*>& lhs, const pair<string, SNP*>& rhs)
{
    return lhs.second->col < rhs.second->col;
}

class Gap {
public:
    uint start;
    uint end;

    Gap(uint s, uint e) {
        start = s;
        end   = e;
    }
};

class Aln {
public:
    uint   id;
    uint   gap_cnt;
    double pct_id;
    string cigar;
    Aln() {
        this->id      = 0;
        this->pct_id  = 0.0;
        this->gap_cnt = 0;
    }
    Aln(uint id, string cigar, double pct_id, uint gap_cnt) {
        this->id      = id;
        this->cigar   = cigar;
        this->pct_id  = pct_id;
        this->gap_cnt = gap_cnt;
    }
};

class PStack {
 public:
    uint            id;
    uint         count; // Number of identical reads forming this stack
    DNANSeq       *seq; // Sequence read
    vector<char *> map; // List of sequence read IDs merged into this stack
    PhyLoc         loc; // Physical genome location of this stack.

    PStack()  {
        id     = 0;
        count  = 0;
        seq    = NULL;
    }
    PStack(const PStack& other);
    PStack& operator= (PStack&& other);
    PStack& operator= (const PStack& other) = delete;

    ~PStack() {
        if (seq!=NULL)
            delete seq;
        for (unsigned int i = 0; i < this->map.size(); i++)
            delete [] this->map[i];
    }

    int  add_id(const char *);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
    void add_read(const char* read_name) {
        char* copy = new char[strlen(read_name)+1];
        strcpy(copy, read_name);
        map.push_back(copy);
        ++count;
    }

    // extend(): Extends the PStack to the desired span.
    void extend(const PhyLoc& phyloc, int length);

    void clear();
    bool operator< (const PStack& other) const;

    static void set_id_of(set<PStack>::iterator pstack, int id) {
        const_cast<PStack&>(*pstack).id = id;
    }
    static void add_read_to(set<PStack>::iterator pstack, const char* read_name) {
        const_cast<PStack&>(*pstack).add_read(read_name);
    }

};

class Stack {
 public:
    uint         id;
    DNANSeq     *seq; // Sequence read
    vector<uint> map; // List of sequence read IDs merged into this stack

    Stack()  {
        id  = 0;
        seq = NULL;
    }
    ~Stack() {
        delete this->seq;
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
};

class Rem {
 public:
    uint          id;
    vector<uint> map; // List of sequence read IDs merged into this stack
    DNANSeq     *seq; // Sequence read
    bool    utilized;

    Rem();
    Rem(int, uint, DNANSeq *);
    ~Rem() {
        delete this->seq;
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
};

class CatMatch {
public:
    int    batch_id;
    int    cat_id; // c-locus ID
    int    sample_id; // sample ID
    int    tag_id; // s-locus ID
    int    depth;
    double lnl;
    char  *haplotype;
    char  *cigar;

    CatMatch() {
        batch_id  = 0;
        cat_id    = 0;
        sample_id = 0;
        tag_id    = 0;
        depth     = 0;
        lnl       = 0.0;
        haplotype = NULL;
        cigar     = NULL;
    }
    ~CatMatch() {
        delete [] haplotype;
        delete [] cigar;
    }
};

class ModRes {
public:
    int   sample_id;
    int   tag_id;
    char *model;

    ModRes(int samp_id, int tag_id, const char *model) {
        this->sample_id = samp_id;
        this->tag_id    = tag_id;
        this->model     = new char[strlen(model) + 1];
        strcpy(this->model, model);
    }
    ~ModRes() {
        delete [] this->model;
    }
};

class SNPRes {
public:
    int   sample_id;
    int   tag_id;
    vector<SNP *> snps;

    SNPRes(int samp_id, int tag_id) {
        this->sample_id = samp_id;
        this->tag_id    = tag_id;
    }
    ~SNPRes() {
        for (uint i = 0; i < this->snps.size(); i++)
            delete this->snps[i];
        this->snps.clear();
    }
};

struct Read {
    DNASeq4 seq;
    string name;

    Read() = default;
    Read(DNASeq4&& s, string&& n)
        : seq(move(s)), name(move(n))
        {}
    Read(Read&&) = default;
    Read& operator= (Read&&) = default;

    bool is_read2() const;
};

struct SiteCounts {
    Counts<Nt2> tot; // The sum over all samples.
    vector<Counts<Nt2>> samples; // With size() == mpopi->samples().size()/
    const MetaPopInfo* mpopi;
};

//
// Inline definitions
// ----------
//

inline
PStack::PStack(const PStack& other)
        : id(other.id)
        , count (other.count)
        , seq (new DNANSeq(*other.seq))
        , map ()
        , loc (other.loc)
        {
    map.reserve(other.map.size());
    for (const char* readname : other.map) {
        char* copy = new char[strlen(readname)+1];
        strcpy(copy, readname);
        map.push_back(copy);
    }
}

inline
PStack& PStack::operator= (PStack&& other) {
    id = other.id;
    count = other.count;
    seq = other.seq;
    other.seq = NULL;
    swap(map, other.map);
    swap(loc, other.loc);
    return *this;
}

inline
void PStack::clear() {
    id = 0;
    count = 0;
    if(seq!=NULL) {
        delete seq;
        seq = NULL;
    }
    map.clear();
    loc.clear();
}

inline
bool PStack::operator< (const PStack& other) const {
    if (loc < other.loc)
        return true;
    else if (other.loc < loc)
        return false;
    else
        // Same genomic loci, compare sequences.
        return *seq < *other.seq;
}

inline
bool Read::is_read2() const {
    return name.length() >= 2
            && name.back() == '2'
            && name[name.length()-2] == '/';
}

#endif // __STACKS_H__
