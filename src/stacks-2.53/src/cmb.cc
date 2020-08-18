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
// cmb.cc -- routines to implement the Combination generating class: CombSet.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: cmb.cc 1990 2010-11-03 04:49:22Z catchen $
//
#include "cmb.h"

using namespace mst;

//
// A cache to store combinations generated as (N choose k)
// for all set sizes encountered
//
map<int, map<int, int **> > _cmbs;

CombSet::~CombSet() {
    int   num_sets, set, i;
    Cmb **c;

    num_sets = this->compound_comb.size();

    for (set = 0; set < num_sets; set++) {
        c = this->compound_comb[set];
        this->destroy(c);
    }

    //
    // Delete edge map
    //
    for (i = 0; i < (int) this->node_list.size(); i++)
        delete [] this->edges[i];
    delete [] this->edges;
}

CombSet::CombSet(int n, int k, MinSpanTree *tree) {
    this->max_set_size = n >= k ? k : n;
    this->num_elements = n - 1;
    this->index        = 0;
    this->mst          = tree;

    int  set_size = this->max_set_size - 1;
    int      size;
    int    **comb;

    cerr << "  Generating combinations for a set of " << n << " elements, with a maximum subset size of " << k << "\n";

    //
    // Add the initial combination: the empty set
    //
    if (_cmbs.count(this->num_elements) == 0 &&
        _cmbs[this->num_elements].count(0) == 0) {
        //cerr << "    N: " << this->num_elements << "; K: 0; Total elements: 0\n";
        comb       = new int * [2];
        comb[0]    = new int[1];
        comb[0][0] = -1;
        comb[1]    = NULL;
        _cmbs[this->num_elements][0] = comb;
    }

    while (set_size > 0) {
        //
        // Check if this set of combinations is already cached.
        //
        if (_cmbs.count(this->num_elements) > 0 &&
            _cmbs[this->num_elements].count(set_size) > 0) {
            set_size--;
            continue;
        }

        //
        // How many combinations will we make?
        //
        size = (int) this->num_combinations(this->num_elements, set_size);

        cerr << "    N: " << this->num_elements << "; K: " << set_size << "; Total elements: " << size << "\n";

        //
        // Generate all combinations, N choose K; N=num_elements, K=set_size
        //
        comb = this->generate_combinations(this->num_elements, set_size, size);

        //
        // Cache this set of combinations
        //
        _cmbs[this->num_elements][set_size] = comb;

        set_size--;
    }

    this->catalog_tree();

    //
    // Finally, generate all combinations of nodes in the tree.
    //
    int max = this->max_set_size < this->num_elements ? this->max_set_size : this->num_elements + 1;
    for (set_size = 0; set_size < max; set_size++)
        this->partition_tree(set_size);

    cerr << "  Total compound combinations for sets of size " << n << ": " << this->compound_comb.size() << "\n";
}

int CombSet::catalog_tree() {
    set<int>      visited;
    queue<Node *> q;
    uint          i, n_1, n_2;

    //
    // Create a two-dimensional array to represent edges between nodes in the tree
    //
    uint cnt   = this->mst->node_count();
    //cerr << "Creating a two-dimensional array of size: " << cnt << " x " << cnt << "\n";
    this->edges = new int * [cnt];
    for (i = 0; i < cnt; i++)
        this->edges[i] = new int[cnt];

    Node *n = this->mst->head();
    q.push(n);
    cnt = 0;

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        this->node_list.push_back(n);
        this->node_map[n->id] = cnt;
        cnt++;

        for (i = 0; i < n->min_adj_list.size(); i++)
            if (visited.count(n->min_adj_list[i]->id) == false)
                q.push(n->min_adj_list[i]);
    }

    n = this->mst->head();
    q.push(n);
    visited.clear();

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (i = 0; i < n->min_adj_list.size(); i++) {

            if (visited.count(n->min_adj_list[i]->id) == false) {
                n_1 = this->node_map[n->id];
                n_2 = this->node_map[n->min_adj_list[i]->id];

                // Create a list of min spanning edges
                this->edge_list.push_back(make_pair(n->id, n->min_adj_list[i]->id));
                // Mark the nodes as connected in our edge array
                this->edges[n_1][n_2] = 1;
                this->edges[n_2][n_1] = 1;
                // Queue this node to be visited next
                q.push(n->min_adj_list[i]);
            }
        }
    }

    return 0;
}

int CombSet::partition_tree(uint set_size) {
    uint          i, j;
    set<int>      visited;
    queue<Node *> q;
    Node         *n;
    int           n_1, n_2, node_cnt, cmb_cnt;
    Cmb         **new_comb, *cmb;
    list<Node *>  nlist_work;

    int **comb     = _cmbs[this->num_elements][set_size];
    int  *subgraph = new int[this->node_list.size()];

    //
    // We want to methodically remove every set of branches of set_size size. The
    // subgraphs represent the combinations we want to generate.
    //
    for (i = 0; comb[i] != NULL; ++i) {
        //
        // This compound combination will consist of set_size+1 subgraphs
        //
        new_comb = new Cmb * [set_size + 2];
        new_comb[set_size + 1] = NULL;

        //
        // Initialize working node list.
        //
        nlist_work = this->node_list;

        //
        // Remove edges
        //
        for (j = 0; j < set_size; j++) {
            n_1 = this->edge_list[comb[i][j]].first;
            n_2 = this->edge_list[comb[i][j]].second;
            this->edges[this->node_map[n_1]][this->node_map[n_2]] = 0;
            this->edges[this->node_map[n_2]][this->node_map[n_1]] = 0;
        }

        //
        // Traverse the subgraphs of the tree and record combinations.
        //
        visited.clear();
        cmb_cnt = 0;
        while (nlist_work.size() > 0) {
            node_cnt = 0;
            n = nlist_work.front();
            q.push(n);
            nlist_work.pop_front();
            //subgraph[node_cnt] = n->id;

            while (!q.empty()) {
                n = q.front();
                q.pop();
                visited.insert(n->id);

                subgraph[node_cnt] = n->id;
                node_cnt++;
                nlist_work.remove(n);

                for (j = 0; j < n->min_adj_list.size(); j++) {
                    n_1 = this->node_map[n->id];
                    n_2 = this->node_map[n->min_adj_list[j]->id];

                    if (visited.count(n->min_adj_list[j]->id) == false &&
                        edges[n_1][n_2] == 1) {
                        q.push(n->min_adj_list[j]);
                    }
                }
            }

            //
            // Package up this combination.
            //
            cmb = new Cmb;
            cmb->size = node_cnt;
            cmb->elem = new int[cmb->size];
            for (j = 0; j < cmb->size; j++)
                cmb->elem[j] = subgraph[j];
            new_comb[cmb_cnt] = cmb;
            cmb_cnt++;
        }

        this->compound_comb.push_back(new_comb);

        //
        // Reset the edges.
        //
        for (j = 0; j < set_size; j++) {
            n_1 = this->edge_list[comb[i][j]].first;
            n_2 = this->edge_list[comb[i][j]].second;
            this->edges[this->node_map[n_1]][this->node_map[n_2]] = 1;
            this->edges[this->node_map[n_2]][this->node_map[n_1]] = 1;
        }
    }

    delete [] subgraph;

    return 0;
}

int **CombSet::generate_combinations(int n, int k, int total) {
    int **comb;

    //
    // Generate an int pointer for each combination, terminate the list with
    // a NULL pointer.
    //
    comb = new int * [total + 1];
    for (int i = 0; i < total; i++)
        comb[i] = new int[k];
    comb[total] = NULL;

    //
    // Setup the initial combination
    //
    int comb_num = 0;

    for (int i = 0; i < k; i++)
        comb[comb_num][i] = i;
    comb_num++;

    //
    // Generate each successive combination
    //
    while (comb_num < total) {
        for (int i = 0; i < k; i++)
            comb[comb_num][i] = comb[comb_num - 1][i];
        this->next_combination(comb[comb_num], n, k);
        comb_num++;
    }

    return comb;
}

int CombSet::next_combination(int *comb, int n, int k) {
    int i;

    //
    // The zero'th position has been incremented to its maximal value,
    // it's not possible to further increment values in the set.
    //
    if (comb[0] > n - k)
        return 0;

    //
    // Increment the last position in the set.
    //
    i = k - 1;
    comb[i]++;

    //
    // Check if the last position has reached its maximal possible value,
    // if so, move back one position, and increment it.
    //
    while ((i > 0) && (comb[i] >= n - k + 1 + i)) {
        i--;
        comb[i]++;
    }

    //
    // Move from the position we incremented above back out to the final position
    //
    for (i = i + 1; i < k; i++)
        comb[i] = comb[i - 1] + 1;

    return 1;
}

long int CombSet::num_combinations(int n, int k) {
    //
    // Compute the binomial coefficient using the method of:
    // Y. Manolopoulos, "Binomial coefficient computation: recursion or iteration?",
    // ACM SIGCSE Bulletin, 34(4):65-67, 2002.
    //
    long int r = 1;
    long int s = (k < n - k) ? n - k + 1 : k + 1;

    for (long int i = n; i >= s; i--)
        r = r * i / (n - i + 1);

    return r;
}

//
// Return a variable length array of Cmb objects, terminated by a NULL pointer.
//
Cmb **CombSet::next(int map[]) {

    if (this->index >= (int) this->compound_comb.size())
        return NULL;

//     int  index, i, j, k, n;
//     int  size = this->compound_comb[this->index]->size;
//     int *e    = this->compound_comb[this->index]->elem;

//     Cmb **c = new Cmb * [size + 1];

//     for (i = 0; i < size; i++) {
//         index = e[i];
//         // sets vector index number
//         k = this->compound_set[index].first;
//         // combination number
//         n = this->compound_set[index].second;

//         c[i] = new Cmb;
//         c[i]->size = this->size[k];
//         c[i]->elem = new int[this->size[k]];

//         for (j = 0; j < this->size[k]; j++)
//             c[i]->elem[j] = (map == NULL) ?
//                 this->sets[k][n][j] :
//                 map[this->sets[k][n][j]];
//     }

//     c[size] = NULL;

    Cmb **c = this->compound_comb[this->index];

    this->index++;

    return c;
}

void CombSet::reset() {
    this->index = 0;
}

void CombSet::destroy(Cmb **cmb) {

    for (uint j = 0; cmb[j] != NULL; j++) {
        delete [] cmb[j]->elem;
        delete cmb[j];
    }
    delete [] cmb;
}

void write_cmb(int *comb, int size) {
    stringstream s;
    string t;

    s << "{";

    for (int i = 0; i < size; i++)
        s << comb[i] << ", ";
    t = s.str().substr(0, s.str().length() - 2);
    t += "}";

    cerr << t << "\n";
}
