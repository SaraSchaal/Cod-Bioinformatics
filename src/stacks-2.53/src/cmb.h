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

#ifndef __CMB_H__
#define __CMB_H__

#include <cmath>
#include <iostream>
#include <sstream>
#include <list>
using std::list;
#include <vector>
#include <string>
#include <set>
#include <utility>

#include "constants.h"
#include "mst.h"
#include "utils.h"

typedef unsigned int uint;

void write_cmb(int *, int);

typedef struct cmb {
    uint size;
    int *elem;
} Cmb;

class CombSet {
    //
    // Given these two variables, we will select N choose K combinations.
    // This combination will be stored in sets, and we will then decrement K by 1
    // and continue to generate sets.
    //
    // Once we have generated all the combinations of a particular size, K, we
    // will partition the minimum spanning tree by dropping combinations of edges
    // from the graph. The selection of edges to drop is provided by the combinations
    // generated first. Finally, each set of disconnected subgraphs makes for one
    // possible combination.
    //
    int num_elements;  // N elements from which we wish to produce combinations
    int max_set_size;  // maximum set size, K, the largest subset we wish to select.

    map<int, int>           node_map;  // Convert non-contiguous IDs from the MST into array indexes for this->edges
    list<mst::Node *>            node_list;
    vector<pair<int, int> > edge_list;
    int                   **edges;

    int            index;
    vector<Cmb **> compound_comb;
    mst::MinSpanTree   *mst;

    int      catalog_tree();
    int      partition_tree(uint);
    int    **generate_combinations(int, int, int);
    int      next_combination(int *, int, int);
    long int num_combinations(int, int);
    void     destroy(Cmb **);

 public:
    CombSet(int, int, mst::MinSpanTree *);
    ~CombSet();

    Cmb **next(int map[] = NULL);
    void  reset();
};

#endif // __CMB_H__
