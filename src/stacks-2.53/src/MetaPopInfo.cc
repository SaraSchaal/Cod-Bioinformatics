#include <fstream>
#include <dirent.h>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>

#include "config.h"
#include "constants.h"
#include "input.h"

#include "MetaPopInfo.h"

using namespace std;

const string Pop::default_name = "defaultpop";
const string Group::default_name = "defaultgrp";

void MetaPopInfo::reset_sample_map() {
    sample_indexes_.clear();
    for (size_t i = 0; i < samples_.size(); ++i) {
        bool rv = sample_indexes_.insert( {samples_[i].name, i} ).second;
        if (!rv) {
            cerr << "Error: Two or more samples have the same name '" << samples_[i].name << "'.\n";
            throw exception();
        }
    }
}

void MetaPopInfo::reset_pop_map() {
    pop_indexes_.clear();
    for (size_t i = 0; i < pops_.size(); ++i) {
        bool rv = pop_indexes_.insert( {pops_[i].name, i} ).second;
        if (!rv) {
            cerr << "Error: Two or more populations have the same name '" << pops_[i].name << "'.\n";
            throw exception();
        }
    }
}

void MetaPopInfo::reset_group_map() {
    group_indexes_.clear();
    for (size_t i = 0; i < groups_.size(); ++i) {
        bool rv = group_indexes_.insert( {groups_[i].name, i} ).second;
        if (!rv) {
            cerr << "Error: Two or more groups have the same name '" << groups_[i].name << "'.\n";
            throw exception();
        }
    }
}

void MetaPopInfo::reset_orig_order() {
    sample_indexes_orig_order_.clear();
    for (const string& name : orig_sample_order_) {
        auto itr = sample_indexes_.find(name);
        if (itr != sample_indexes_.end())
            // This sample is still there (n.b. delete_samples/intersect_with).
            sample_indexes_orig_order_.push_back(itr->second);
    }
    assert(sample_indexes_orig_order_.size() == samples_.size());
}

void MetaPopInfo::init_popmap(const string& pmap_path) {
    assert(samples_.empty());

    ifstream fh(pmap_path.c_str(), ifstream::in);
    if (fh.fail()) {
        cerr << "Error: Failed to open population map file '" << pmap_path << "'.\n";
        throw exception();
    }

    size_t p = 0; // pop index counter
    size_t g = 0; // group index counter

    char line[max_len];
    memset(line, '\0', max_len);
    vector<string> parts;
    while (fh.getline(line, max_len)) {
        size_t len = strlen(line);

        // Skip empty lines and comments
        if (len == 0 || line[0] == '#')
            continue;

        // Check for Windows line endings
        if (line[len - 1] == '\r') {
            line[len - 1] = '\0';
            len -= 1;
        }

        //
        // Parse the contents, we expect:
        // <file name> tab <population string> [tab <group string>]
        //

        parse_tsv(line, parts);
        for (string& part : parts)
            strip_whitespace(part);

        if (parts.size() < 2
            || parts.size() > 3
            || parts[0].empty()
            || parts[1].empty()) {
            cerr << "Error: Malformed population map -- expected 'SAMPLE \\t POP [\\t GROUP]'. In file '" << pmap_path << "', line :\n" << line << "\n";
            throw exception();
        }

        //
        // Process the sample name.
        //

        samples_.push_back(Sample(parts[0]));
        orig_sample_order_.push_back(parts[0]);

        //
        // Process the population field.
        //

        pair<map<string,size_t>::iterator, bool> pop_ins = pop_indexes_.insert( {parts[1], p} );
        size_t pop_index = pop_ins.first->second;

        samples_.back().pop = pop_index; // Set the sample's population index.
        if (pop_ins.second) {
            // Unknown pop
            pops_.push_back(Pop(parts[1]));
            ++p;
        }

        //
        // Process the group field, if any
        //

        if (parts.size() == 3 && ! parts[2].empty()) {
            // Get the index of this group -- create it if necessary.
            pair<map<string,size_t>::iterator, bool> grp_ins = group_indexes_.insert( {parts[2], g} );
            if (grp_ins.second) {
                groups_.push_back(Group(parts[2]));
                ++g;
            }
            size_t grp_index = grp_ins.first->second;

            if (pops_[pop_index].group != size_t(-1)) {
                // The current pop already has a group, check that it is the same.
                if (pops_[pop_index].group != grp_index) {
                    cerr << "Error: In population map file '"
                         << pmap_path << "': population '"
                         << pops_[pop_index].name << "' belongs to two groups, '"
                         << groups_[pops_[pop_index].group].name << "' and '"
                         << groups_[grp_index].name << "'\n.";
                    throw exception();
                }
            } else {
                pops_[pop_index].group = grp_index;
                groups_[grp_index].pops.push_back(pop_index);
            }
        }
    }
    if (samples_.empty()) {
        cerr << "Error: Population map '" << pmap_path << "' appears to be empty.\n";
        throw exception();
    }

    //
    // Check that all the populations are in a group. Put
    // populations that do not have a group in a default
    // one.
    //

    bool missing_group = false;
    if (not groups_.empty()) {
        for (vector<Pop>::iterator p = pops_.begin(); p != pops_.end(); ++p) {
            if (p->group == size_t(-1)) {
                cerr << "Warning: Population '" << p->name
                << "' did not have a group, adding it to '"
                << Group::default_name << "'.\n";
                missing_group = true;
            }
        }
    } else {
        missing_group = true;
    }
    if (missing_group) {
        groups_.push_back(Group(Group::default_name));
        g = groups_.size()-1;
        group_indexes_.insert( {Group::default_name, g} );
        for (size_t p = 0; p < pops_.size(); ++p) {
            if (pops_[p].group == size_t(-1)) {
                pops_[p].group = g;
                groups_[g].pops.push_back(p);
            }
        }
    }

    //
    // Sort the samples. Determine the first/last indexes for each population.
    //

    sort(samples_.begin(), samples_.end());
    reset_sample_map();
    reset_orig_order();

    size_t curr_pop = 0;
    pops_[curr_pop].first_sample = 0;
    for (size_t s = 1; s < samples_.size(); ++s) {
        if (samples_[s].pop != curr_pop) {
            pops_[curr_pop].last_sample = s-1;
            ++curr_pop;
            pops_[curr_pop].first_sample = s;
        }
    }
    pops_[curr_pop].last_sample = samples_.size()-1;
}

void MetaPopInfo::init_names(const vector<string>& sample_names) {
    if(sample_names.empty()) {
        cerr << "Error: No samples.\n";
        throw exception();
    }
    assert(samples_.empty());

    orig_sample_order_ = sample_names;

    // Create the samples
    for (vector<string>::const_iterator s=sample_names.begin(); s!= sample_names.end(); ++s) {
        samples_.push_back(Sample(*s));
        samples_.back().pop = 0;
    }
    sort(samples_.begin(), samples_.end());

    // Create a default population
    pops_.push_back(Pop(Pop::default_name));
    pops_[0].first_sample = 0;
    pops_[0].last_sample = samples_.size()-1;
    pops_[0].group = 0;

    // Create a default group
    groups_.push_back(Group(Group::default_name));
    groups_[0].pops.push_back(0);

    // Set the support members.
    reset_sample_map();
    reset_pop_map();
    reset_group_map();
    reset_orig_order();
}

void MetaPopInfo::init_directory(const string& dir_path) {

    //
    // Find all sample names.
    //

    vector<string> sample_names;
    DIR* dir = opendir(dir_path.c_str());
    if (dir == NULL) {
        cerr << "Unable to open directory '" << dir_path << "' for reading.\n";
        throw exception();
    }
    dirent *direntry;
    while ((direntry = readdir(dir)) != NULL) {
        string filename = direntry->d_name;

        if (filename == "." || filename == ".." || filename.substr(0, 6) == "batch_")
            continue;

        size_t pos = filename.rfind(".tags.tsv");
        if (pos == string::npos)
            pos = filename.rfind(".tags.tsv.gz");

        if (pos != string::npos)
            sample_names.push_back(filename.substr(0, pos));
    }
    closedir(dir);

    if (sample_names.empty()) {
        cerr << "Error: Failed to find sample files in directory '" << dir_path << "'.\n";
        throw exception();
    }

    //
    // Initialize the MetaPopInfo
    //
    init_names(sample_names);
}

void MetaPopInfo::delete_samples(const vector<size_t>& rm_samples) {

    if (rm_samples.empty())
        return;

    // Remove these samples from [samples_].
    for (size_t s : rm_samples)
        samples_.at(s).name.clear(); // Mark the sample for removal.
    samples_.erase(remove_if(
            samples_.begin(),
            samples_.end(),
            [] (Sample& s) {return s.name.empty();}
            ), samples_.end());

    // Update the indexes of the populations.
    for (Pop& p : pops_) {
        for (vector<size_t>::const_reverse_iterator rm_sample = rm_samples.rbegin(); rm_sample != rm_samples.rend(); ++rm_sample) {
            if (p.first_sample > *rm_sample) // n.b. ">"
                --p.first_sample;
            if (p.last_sample >= *rm_sample) // n.b. ">=". Thus if the population becomes
                                              // empty, [first_sample] will be past [last_sample].
                                              // n.b. If removing the first pop, last_sample=size_t(-1).
                --p.last_sample;
        }
    }

    // Remove empty populations.
    auto pop_is_empty = [] (Pop& p) {return (p.first_sample > p.last_sample || p.last_sample == size_t(-1));};
    for(Group& g : groups_)
        g.pops.erase(remove_if(
                g.pops.begin(),
                g.pops.end(),
                [this,&pop_is_empty] (size_t p) {return pop_is_empty(pops_[p]);}
                ), g.pops.end());
    pops_.erase(remove_if(
            pops_.begin(),
            pops_.end(),
            pop_is_empty
            ),  pops_.end());

    // Remove empty groups from [groups_].
    groups_.erase(remove_if(
            groups_.begin(),
            groups_.end(),
            [] (Group& g) {return g.pops.empty();}
            ), groups_.end());

    // Update the support members.
    reset_sample_map();
    reset_pop_map();
    reset_group_map();
    reset_sample_id_map();
    reset_orig_order();
}

void MetaPopInfo::intersect_with(const vector<string>& samples) {
    vector<size_t> common_samples;
    for (const string& s : samples) {
        auto itr = sample_indexes_.find(s);
        if (itr != sample_indexes_.end())
            common_samples.push_back(itr->second);
    }
    sort(common_samples.begin(), common_samples.end());

    vector<size_t> rm_samples;
    auto next_common = common_samples.begin();
    for (size_t i=0; i< samples_.size(); ++i)
        if (next_common != common_samples.end() && i == *next_common)
            ++next_common;
        else
            rm_samples.push_back(i);

    delete_samples(rm_samples);
}

void
MetaPopInfo::status(ostream &fh)
{
    fh << "Working on " << this->samples().size() << " samples.\n";
    fh << "Working on " << this->pops().size() << " population(s):\n";
    for (vector<Pop>::const_iterator p = this->pops().begin(); p != this->pops().end(); p++) {
        size_t indent   = p->name.length() + 6;
        size_t line_lim = 120;
        size_t line_len = 0;

        fh << "    " << p->name << ": ";
        for (size_t s = p->first_sample; s < p->last_sample; ++s) {
            fh << this->samples()[s].name << ", ";
            line_len += this->samples()[s].name.length() + 2;
            if (line_len > line_lim) {
                fh << "\n" << string(indent, ' ');
                line_len = 0;
            }
        }
        fh << this->samples()[p->last_sample].name << "\n";
    }
    fh << "Working on " << this->groups().size() << " group(s) of populations:\n";
    for (vector<Group>::const_iterator g = this->groups().begin(); g != this->groups().end(); g++) {
        fh << "    " << g->name << ": ";
        for (vector<size_t>::const_iterator p = g->pops.begin(); p != g->pops.end() -1; ++p) {
            //rem. end()-1 and back() are safe, there's always at least one pop
            fh << this->pops()[*p].name << ", ";
        }
        fh << this->pops()[g->pops.back()].name << "\n";
    }
    fh << "\n";
}

void MetaPopInfo::reset_sample_id_map() {
    sample_indexes_by_id_.clear();
    for (size_t i = 0; i < samples_.size(); ++i)
        sample_indexes_by_id_[samples_[i].id] = i;
}

void MetaPopInfo::fill_files(vector<pair<int, string> >& files) const {
    files.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        files.push_back( {sample->pop, sample->name} );
}

void MetaPopInfo::fill_sample_ids(vector<int>& sample_ids) const {
    sample_ids.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        sample_ids.push_back(sample->id);
}

void MetaPopInfo::fill_samples(map<int, string>& samples) const {
    samples.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        samples.insert({sample->id, sample->name});
}

void MetaPopInfo::fill_pop_key(map<int, string>& pop_key) const {
    pop_key.clear();
    for (size_t i = 0; i < pops_.size(); ++i)
        pop_key.insert( {i, pops_[i].name} );
}

void MetaPopInfo::fill_pop_indexes(map<int, pair<int, int> >& pop_indexes) const {
    pop_indexes.clear();
    for (size_t i = 0; i < pops_.size(); ++i)
        pop_indexes.insert( {i, {pops_[i].first_sample, pops_[i].last_sample}} );
}

void MetaPopInfo::fill_grp_key(map<int, string>& grp_key) const {
    grp_key.clear();
    for (size_t i = 0; i < groups_.size(); ++i)
        grp_key.insert({i, groups_[i].name});
}

void MetaPopInfo::fill_grp_members(map<int, vector<int> >& grp_members) const {
    grp_members.clear();
    for (size_t i = 0; i < groups_.size(); ++i) {
        vector<int>& pop_ids = grp_members.insert( {i, vector<int>()} ).first->second;
        for(vector<size_t>::const_iterator p = groups_[i].pops.begin(); p != groups_[i].pops.end(); ++p)
            pop_ids.push_back(*p);
    }
}
