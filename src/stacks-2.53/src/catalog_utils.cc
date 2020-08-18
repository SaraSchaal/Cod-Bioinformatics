// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
// catalog_utils.cc -- common routines for manipulating catalog objects.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include <regex>

#include "utils.h"

#include "catalog_utils.h"

using namespace std;

const regex catalog_tags_regex ("^batch_([0-9]+).catalog.tags.tsv(.gz)?$");

vector<int> find_catalogs(const string& dir_path) {
    vector<int> ids;

    for (DirIterator e (dir_path); e; ++e) {
        smatch m;
        string name (e.name());
        regex_match(name, m, catalog_tags_regex);
        if (!m.empty())
            ids.push_back(stoi(m[1].str()));
    }

    return ids;
}

int
reduce_catalog(map<int, CSLocus *> &catalog, set<int> &whitelist, set<int> &blacklist)
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0)
        return 0;

    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
        if (blacklist.count(loc->id)) continue;

        list[it->first] = it->second;
        i++;
    }

    catalog = list;

    return i;
}

int
check_whitelist_integrity(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist)
{
    if (whitelist.size() == 0) return 0;

    int rm_snps = 0;
    int rm_loci = 0;

    CSLocus *loc;
    map<int, set<int> >::iterator it;
    set<int>::iterator sit;
    map<int, set<int> > new_wl;

    cerr << "Checking the integrity of the whitelist...";

    for (it = whitelist.begin(); it != whitelist.end(); it++) {
        if (catalog.count(it->first) == 0) {
            rm_loci++;
            cerr << "\n  Removing locus " << it->first << " from whitelist as it does not exist in the catalog.";
        } else {
            loc = catalog[it->first];

            if (it->second.size() == 0) {
                new_wl.insert(make_pair(it->first, set<int>()));
                continue;
            }

            set<int> cat_snps;
            for (uint i = 0; i < loc->snps.size(); i++)
                cat_snps.insert(loc->snps[i]->col);

            for (sit = it->second.begin(); sit != it->second.end(); sit++)
                if (cat_snps.count(*sit)) {
                    new_wl[it->first].insert(*sit);
                } else {
                    rm_snps++;
                    cerr << "\n  Removing SNP at column " << *sit << " in locus " << it->first << " from whitelist as it does not exist in the catalog.";
                }
        }
    }

    whitelist = new_wl;

    if (rm_loci > 0 || rm_snps > 0) cerr << "\n";

    cerr << "done.\n"
         << "Removed " << rm_loci << " loci and " << rm_snps << " SNPs from the whitelist that were not found in the catalog.\n";

    return 0;
}

int
reduce_catalog(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist, set<int> &blacklist)
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0)
        return 0;

    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
        if (blacklist.count(loc->id)) continue;

        list[it->first] = it->second;
        i++;
    }

    catalog = list;

    return i;
}
