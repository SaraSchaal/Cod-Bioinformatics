// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2018, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __MSTACK_H__
#define __MSTACK_H__

#include <string>
#include <vector>
#include <map>
#include <utility>
#include<iostream>

#include "stacks.h"

class MergedStack {
 public:
    int   id;     // Identifier for the merged stack.
    char *con;    // Consensus sequence
    uint  len;    // Sequence length

    //
    // Stack component parts
    //
    int                       count; // Number of merged stacks
    vector<int>               utags; // Stack IDs that have been merged into this MergedStack
    vector<int>             remtags; // Remainder tag IDs that have been merged into this Stack
    DNANSeq                **matrix; // Two-dimensional array for iterating over the combined stack (stacks and remainders).
    vector<pair<int, int> >    dist; // Vector describing the distance between this stack and other stacks.
    vector<Aln>                alns; // Vector describing gapped alignments between this stack and other stacks.
    vector<size_t>        rem_queue; // Temporary list of Remainder reads that are queue for merging.
    
    int cohort_id; // Group ID of all stacks that were originally part of the same subgraph
    double    lnl; // Log likelihood of this stack

    //
    // Mapping components
    //
    PhyLoc               loc; // Physical genome location of this Stack.
    vector<SNP *>       snps; // Single Nucleotide Polymorphisms found in this Stack
    map<string, int> alleles; // Set of alleles defined by the SNPs found in this Stack
    vector<Gap>         gaps;

    //
    // Flags
    //
    bool deleveraged;
    bool masked;
    bool blacklisted;
    bool gappedlumberjack;
    bool lumberjackstack;

    MergedStack();
    ~MergedStack();

    MergedStack(MergedStack&& other);
    MergedStack(const MergedStack& other) = delete;
    MergedStack& operator=(const MergedStack& other) = delete;

    int       add_consensus(const char *);
    int       add_consensus(DNASeq *);
    int       add_consensus(const DNANSeq *);
    int       add_dist(const int id, const int dist);
    DNANSeq **gen_matrix(map<int, Stack *> &, map<int, Rem *> &);
    DNANSeq **gen_matrix(map<int, PStack *> &);
    double    calc_likelihood();
    string    write_cmb();
};

inline
MergedStack::MergedStack(MergedStack&& o)
        : id (o.id)
        , con (o.con)
        , len (o.len)
        , count (o.count)
        , utags (move(o.utags))
        , remtags (move(o.remtags))
        , matrix (o.matrix)
        , dist (move(o.dist))
        , alns (move(o.alns))
        , cohort_id (o.cohort_id)
        , lnl (o.lnl)
        , loc (move(o.loc))
        , snps (move(o.snps))
        , alleles (move(o.alleles))
        , gaps (move(o.gaps))
        , deleveraged (o.deleveraged)
        , masked (o.masked)
        , blacklisted (o.blacklisted)
        , gappedlumberjack (o.gappedlumberjack)
        , lumberjackstack (o.lumberjackstack)
        {
    o.con = NULL;
    o.matrix = NULL;
    o.snps.clear();
}
#endif // __MSTACK_H__
