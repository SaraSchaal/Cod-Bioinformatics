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

#ifndef __MST_H__
#define __MST_H__

#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <queue>
using std::queue;
#include <sstream>
#include<iostream>

#include "constants.h"

namespace mst {

class Node;

class Edge {
 public:
    uint  dist;   // Distance or weight
    Node *child;
};

class Node {
public:
    uint           id;
    string         label;
    vector<Edge *> edges;

    Node *parent;
    bool  update;
    uint  min_dist;

    //
    // List of adjacent nodes that are connected by minimal distance
    //
    vector<Node *> min_adj_list;

    Node(uint id) {
        this->id       = id;
        this->parent   = NULL;
        this->update   = true;
        this->min_dist = 1000000;
    }

    ~Node()  {
        this->min_adj_list.clear();

        for (uint i = 0; i < this->edges.size(); i++)
            delete this->edges[i];
        this->edges.clear();
    }

    Edge *add_edge(Node *, int);
};

class MinSpanTree {
    map<int, Node *> nodes;
    map<string, int> node_key;
    uint             id_cnt;

 public:
    MinSpanTree()  { id_cnt = 0; }
    ~MinSpanTree() {
        for (map<int, Node *>::iterator i = this->nodes.begin(); i != this->nodes.end(); i++)
            delete i->second;
        this->nodes.clear();
    }

    Node  *add_node(int id);
    Node  *add_node(string label);
    int    build_tree();
    int    node_count();
    Node  *node(int id);
    Node  *node(string label);
    Node  *head();
    bool   connected(int *, int);
    string vis(bool);
};

bool min_span_tree_cmp(const Node *, const Node *);

} // namespace mst
#endif // __MST_H__
