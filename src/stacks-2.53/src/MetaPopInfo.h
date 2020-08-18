#ifndef METAPOPINFO_H
#define METAPOPINFO_H

#include <string>
#include <vector>
#include <map>

#include "constants.h"

struct Sample {
    string name;
    size_t pop;
    int id; // optional

    Sample(const string& n) : name(n), pop(SIZE_MAX), id(-1) {}
    inline bool operator<(const Sample& other) const;
};

struct Pop {
    string name;
    size_t first_sample;
    size_t last_sample;
    size_t group;

    Pop(const string& n) : name(n), first_sample(SIZE_MAX), last_sample(SIZE_MAX), group(SIZE_MAX) {}
    size_t n_samples() const {return last_sample - first_sample + 1;}
    static const string default_name;
};

struct Group {
    string name;
    vector<size_t> pops;

    Group(const string& n) : name(n), pops() {}
    static const string default_name;
};

/*
 * MetaPopInfo
 * Class for reprensenting a metapopulation : its samples, populations,
 * groups of populations, and associated information.
 */
class MetaPopInfo {
    vector<Sample> samples_; //n.b. Samples are sorted primarily by population index, and secondarily by name.
    vector<Pop> pops_;
    vector<Group> groups_;

    map<string,size_t> sample_indexes_; // Links a name with an index in [samples_].
    map<string,size_t> pop_indexes_;
    map<string,size_t> group_indexes_;
    void reset_sample_map(); // Resets [sample_indexes_].
    void reset_pop_map();
    void reset_group_map();

    vector<string> orig_sample_order_; // Safeguards the original input sample order.
    vector<size_t> sample_indexes_orig_order_;
    void reset_orig_order();

    map<size_t,size_t> sample_indexes_by_id_; // Links a sample ID with an index in [samples_].
    void reset_sample_id_map();

    MetaPopInfo(MetaPopInfo&& other) = delete; // Immovable (for pointer stability).
public:
    MetaPopInfo() = default;

    // Create the representation :
    // -- from a population map file.
    // -- from just a vector of sample names.
    // -- or by looking for "*.tags.tsv(.gz)" files in a directory.
    void init_popmap(const string& popmap_path);
    void init_names(const vector<string>& sample_names);
    void init_directory(const string& dir_path);

    // Delete samples from the metapopulation.
    // (As samples, populations or groups may be deleted, the indexes of
    // the remaining ones change, but the order in which they appear
    // is preserved.)
    void delete_samples(const vector<size_t>& rm_samples);

    // Intersects the population map with a list of samples.
    // May call `delete_samples()`.
    void intersect_with(const vector<string>& samples);

    // Retrieve information.
    const vector<Sample>& samples() const {return samples_;}
    size_t n_samples() const {return samples().size();}
    const vector<Pop>& pops() const {return pops_;}
    const vector<Group>& groups() const {return groups_;}

    const vector<size_t>& sample_indexes_orig_order() const {return sample_indexes_orig_order_;}

    size_t get_sample_index(const string& name, bool must_exist=true) const;
    size_t get_pop_index(const string& name) const {return pop_indexes_.at(name);}
    size_t get_group_index(const string& name) const {return group_indexes_.at(name);}

    // Work with sample IDs. (IDs unicity is not enforced.)
    void set_sample_id(size_t index, size_t id) {samples_.at(index).id = id; sample_indexes_by_id_[id] = index;}
    size_t get_sample_index(const size_t& id) const {return sample_indexes_by_id_.at(id);}

    void status(ostream &fh);

    /*
     * Methods for backwards compatibility
     */

    // Fill former globals.
    void fill_files(vector<pair<int, string> >&) const;
    void fill_sample_ids(vector<int>&) const;
    void fill_samples(map<int, string>&) const;
    void fill_pop_key(map<int, string>&) const;
    void fill_pop_indexes(map<int, pair<int, int> >&) const;
    void fill_grp_key(map<int, string>&) const;
    void fill_grp_members(map<int, vector<int> >&) const;
};

inline
bool Sample::operator<(const Sample& other) const {
    if (pop == other.pop)
        return name < other.name;
    else
        return pop < other.pop;
}

inline
size_t MetaPopInfo::get_sample_index(
    const string& name,
    bool must_exist
) const {
    if (must_exist) {
        return sample_indexes_.at(name);
    } else {
        auto itr = sample_indexes_.find(name);
        return (itr == sample_indexes_.end() ? SIZE_MAX : itr->second);
    }
}

#endif // METAPOPINFO_H
