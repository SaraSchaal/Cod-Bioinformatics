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
// mst.cc -- routines to implement the Minimum Spanning Tree Class:.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "mst.h"

namespace mst {

Edge *Node::add_edge(Node *n, int dist) {
    Edge *e = new Edge;
    e->child = n;
    e->dist  = dist;
    this->edges.push_back(e);

    return e;
}

bool min_span_tree_cmp(const Node *lhs, const Node *rhs) {
    return (lhs->min_dist > rhs->min_dist);
}

Node *MinSpanTree::add_node(int id) {
    Node *n = new Node(id);

    if (this->nodes.count(id) > 0)
        delete this->nodes[id];

    this->nodes[id] = n;

    return n;
}

Node *MinSpanTree::add_node(string label) {
    //
    // Obtain an ID for this node.
    //
    uint id = this->id_cnt;

    if (this->nodes.count(id) > 0) {
        cerr << "Error constructing minimum spanning tree.\n";
        delete this->nodes[id];
    }

    Node *n = new Node(id);
    n->label = label;

    this->nodes[id]       = n;
    this->node_key[label] = id;
    this->id_cnt++;

    return n;
}

Node *MinSpanTree::node(int id) {
    return this->nodes[id];
}

Node *MinSpanTree::node(string label) {
    uint id = this->node_key[label];
    return this->nodes[id];
}

Node *MinSpanTree::head() {
    return this->nodes.begin()->second;
}

int MinSpanTree::node_count() {
    return this->nodes.size();
}

bool MinSpanTree::connected(int *ids, int size) {
    set<int>      valid, visited;
    queue<Node *> q;

    if (size == 1)
        return true;

    for (int i = 0; i < size; i++)
        valid.insert(ids[i]);

    //
    // Take the first ID and begin traversing the tree. If we hit
    // a node not in the ids set, stop traversing this branch. Check that
    // all nodes are directly connected.
    //
    int   valid_cnt = 0;
    Node *n         = this->node(ids[0]);
    q.push(n);

    while (!q.empty() && valid_cnt < size) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        if (valid.count(n->id)) {
            valid_cnt++;

            for (uint i = 0; i < n->min_adj_list.size(); i++)
                if (visited.count(n->min_adj_list[i]->id) == false)
                    q.push(n->min_adj_list[i]);
        }
    }

    if (valid_cnt == size)
        return true;
    else
        return false;
}

//
// Build a minimum spanning tree using Prim's alogorithm. Assume all necessary
// nodes have been added using the add_node function.
//
int MinSpanTree::build_tree() {
    //
    // Vector, which we treat as a binary heap to access nodes that are of minimal distance
    //
    vector<Node *> q;

    //
    // Select an initial node to process and initialize its minimum distance.
    //
    Node *n = this->nodes.begin()->second;
    n->min_dist = 0;

    //
    // Add all of the nodes to the binary heap; process them in order of min_dist
    //
    map<int, Node *>::iterator it;
    for (it = this->nodes.begin(); it != this->nodes.end(); it++)
        q.push_back((*it).second);
    make_heap(q.begin(), q.end(), min_span_tree_cmp);

    while (q.size() > 0) {
        n = q.front();
        pop_heap(q.begin(), q.end());
        q.pop_back();
        n->update = false;

        //cerr << "Examining node: " << n->id << " (" << n->min_dist << ")\n";

        //
        // Record the minimum connection between parent and n.
        //
        if (n->parent != NULL) {
            n->parent->min_adj_list.push_back(n);
            n->min_adj_list.push_back(n->parent);
        }

        //
        // Iterate through all of the edges of n and update the
        // minimum distance to the proper nodes.
        //
        Edge *e;
        for (uint i = 0; i < n->edges.size(); i++) {
            e = n->edges[i];

            if (e->child->update == true && e->dist < e->child->min_dist) {
                e->child->parent   = n;
                e->child->min_dist = e->dist;

                //cerr << "  Updating node: " << e->child->id << " to have distance: " << e->child->min_dist << "\n";
            }
        }

        //
        // Resort the heap after possibly changing many min_dist values
        //
        make_heap(q.begin(), q.end(), min_span_tree_cmp);
    }

    return 0;
}

string MinSpanTree::vis(bool overlay) {
    uint   j;
    double d, scale, scaled_d;
    char   label[32];
    int    scale_factor = 20;

    //
    // Output a specification to visualize the minimum spanning tree using graphviz:
    //   http://www.graphviz.org/
    //
    stringstream data;
    data << "graph stacks_" << this->nodes.size() << " {\n"
         << "rankdir=LR\n"
         << "size=\"" << scale_factor << "!\"\n"
         << "overlap=false\n"
         << "node [shape=circle style=filled fillcolor=\"#3875d7\" fontname=\"Arial\"];\n"
         << "edge [fontsize=8.0 fontname=\"Arial\" color=\"#aaaaaa\"];\n";

    map<int, Node *>::iterator i;
    set<int>      visited;
    queue<Node *> q;

    //
    // If overlay==true, write the minimum spanning tree on top of the full tree as a subgraph.
    //
    data << "subgraph mst {\n"
         << "    edge [penwidth=5 fontsize=12.0 fontcolor=\"black\" color=\"black\"]\n"
         << "    node [fillcolor=\"red\" fontcolor=\"white\"]\n";

    Node *n = this->head();
    q.push(n);

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            data << "  ";
            n->label.length() > 0 ? data << n->label : data << n->id;
            data << "--";
            n->min_adj_list[i]->label.length() > 0 ? (data << n->min_adj_list[i]->label) : (data << n->min_adj_list[i]->id);
            data << "\n";
            if (visited.count(n->min_adj_list[i]->id) == 0)
                q.push(n->min_adj_list[i]);
        }
    }

    data << "}\n";

    //
    // Scale the graph to display on a scale_factor inch canvas. Find the largest edge weight
    // and scale the edge lengths to fit the canvas.
    //
    scale = 0.0;
    for (i = this->nodes.begin(); i != this->nodes.end(); i++) {
        n = i->second;
        for (j = 0; j < n->edges.size(); j++)
            scale = n->edges[j]->dist > scale ? n->edges[j]->dist : scale;
    }
    scale = scale / scale_factor;

    //
    // Write out edges.
    //
    for (i = this->nodes.begin(); i != this->nodes.end(); i++) {
        n = i->second;
        for (j = 0; j < n->edges.size(); j++) {
            d        = n->edges[j]->dist;
            scaled_d = d / scale;
            scaled_d = scaled_d < 0.75 ? 0.75 : scaled_d;
            sprintf(label, "%.1f", d);

            n->label.length() > 0 ? (data << n->label) : (data << n->id);
            data << " -- ";
            n->edges[j]->child->label.length() > 0 ? (data << n->edges[j]->child->label) : (data << n->edges[j]->child->id);
            data << " [len=" << scaled_d << ", label=" << label << "];\n";
        }
    }

    data << "}\n";

    return data.str();
}

}//namespace
