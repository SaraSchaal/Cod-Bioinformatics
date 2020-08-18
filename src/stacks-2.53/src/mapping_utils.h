// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2019, Julian Catchen <jcatchen@illinois.edu>
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
#ifndef MAPPING_UTILS_H
#define MAPPING_UTILS_H

#include <map>
#include <vector>
#include <string>

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "MetaPopInfo.h"
#include "genotype_dictionaries.h"


//
// MappingGenotypesProcessor
// ----------
// Class for processing genetic markers to derive mappable genotypes.
//
class MappingGenotypesProcessor {
    MetaPopInfo *_mpopi;         // Population Map
    size_t       _progeny_index; // Index in the population map for progeny
    CrossT       _cross_type;    // Type of mapping cross specified

    // Counters
    size_t   _num_loci;
    size_t   _mappable_loci;
    double   _mean_progeny_sum;
    size_t   _corrected;

    map<CrossT, string> _types =
        {{CrossT::cp,  "CP"},
         {CrossT::dh,  "DH"},
         {CrossT::bc1, "BC1"},
         {CrossT::f2,  "F2"}};

    //
    // Translation dictionaries for markers and genotypes.
    //
    map<string, string>              _marker_map;
    map<string, map<string, string>> _genotype_map;
    map<string, map<string, double>> _segregation_ratios;

    //
    // Sets of segregation ratios for detecting and correcting missing alleles in CP maps.
    //
    // For the particular marker type, if we see the following genotype ratios, we
    // can infer that one of the parents is missing an allele.
    //
    map<string, map<string, double>> _c_ratios_1 =  // Correction segregation ratios 1
        {{"ab/aa",
          {{"aa", 0.50},
           {"ab", 0.25},
           {"bb", 0.25}}},
         {"aa/ab",
          {{"aa", 0.50},
           {"ab", 0.25},
           {"bb", 0.25}}},
         {"ab/cc",
          {{"ac", 0.25},
           {"bc", 0.25},
           {"aa", 0.25},
           {"bb", 0.25}}},
         {"cc/ab",
           {{"ac", 0.25},
            {"bc", 0.25},
            {"aa", 0.25},
            {"bb", 0.25}}}};
    map<string, map<string, double>> _c_ratios_2 =  // Correction segregation ratios 2
        {{"aa/bb",
          {{"aa", 0.50},
           {"ab", 0.50}}}};
    map<string, map<string, double>> _c_ratios_3 =  // Correction segregation ratios 3
        {{"aa/bb",
          {{"bb", 0.50},
           {"ab", 0.50}}}};

public:
    MappingGenotypesProcessor(): _mpopi(NULL), _num_loci(0), _mappable_loci(0), _mean_progeny_sum(0.0) {}
    MappingGenotypesProcessor(MetaPopInfo *popi, CrossT ct);

    bool   is_mapping_cross(ostream &log_fh);
    int    next_batch(const vector<LocBin *> &loci);
    int    correct_cp_marker(CSLocus *, Datum **, size_t, size_t, map<string, int> &);
    
    int    type(CrossT ct)  { this->_cross_type = ct; return 0; }
    CrossT type() const     { return this->_cross_type; }
    string type_str() const { return this->_types.at(this->_cross_type); }

    size_t n_loci() const        { return this->_num_loci; }
    size_t mappable_loci() const { return this->_mappable_loci; }
    double mean_progeny() const  { return this->_mean_progeny_sum / this->_mappable_loci; }
    size_t corrected() const     { return this->_corrected; }
};

//
// Chi-squared distribution critical values
//   [0] => p-values
//   [1] => one degree of freedom
//   [2] => two degrees of freedom
//   [3] => three degrees of freedom
//
const int chisq_crit_values_size = 10;
const double chisq_crit_values[4][10] = {
    {0.50, 0.25, 0.20, 0.15, 0.10, 0.05,  0.01,  0.005,  0.001,  0.0005},
    {0.46, 1.32, 1.64, 2.07, 2.71, 3.84,  6.63,  7.880, 10.830, 12.1200},
    {1.39, 2.77, 3.22, 3.79, 4.61, 5.99,  9.21, 10.600, 13.820, 15.2000},
    {2.37, 4.11, 4.64, 5.32, 6.25, 7.81, 11.34, 12.840, 16.270, 17.7300}
};
const double chisq_pval_limit = 0.05;

string classify_marker(const LocBin *, const MetaPopInfo *);
int    assign_parent_genotypes(const MetaPopInfo *, const LocBin *, string, map<string, string> &);
string assign_generic_genotype(string, const Datum *, map<string, string> &);
double chisq_test(map<string, map<string, double>> &, map<string, int> &, string, double);
double chisq_pvalue(int, double);


#endif // MAPPING_UTILS_H
