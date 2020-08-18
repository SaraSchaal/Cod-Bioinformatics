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

#ifndef __PHASEDSTACKS_H__
#define __PHASEDSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>

#include "constants.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "input.h"
#include "sql_utilities.h"
#include "locus.h"

#ifdef HAVE_LIBZ
#include <cerrno>
#include <zlib.h>
#endif

enum loc_t {strong_ld, recomb, uninformative};

class PhasedSample {
public:
    string name;
    int    id;
    int    size;
    char  *nucs_1;
    char  *nucs_2;

    PhasedSample() {
        this->id   = 0;
        this->size = 0;
        this->nucs_1 = NULL;
        this->nucs_2 = NULL;
    }
    ~PhasedSample() {
        if (this->nucs_1 != NULL)
            delete [] this->nucs_1;
        if (this->nucs_2 != NULL)
            delete [] this->nucs_2;
    }
};

class NucSum {
public:
    uint  bp;
    uint  col;
    uint  clocus;
    float freq;
    uint  nuc[4];
    // nuc[0] == A
    // nuc[1] == C
    // nuc[2] == G
    // nuc[3] == T

    NucSum() {
        this->freq    = 1.0;
        this->bp      = 0;
        this->clocus  = 0;
        for (uint i = 0; i < 4; i++)
            this->nuc[i] = 0;
    }
};

class dPrime {
public:
    double dprime;
    bool   chisq_p;
    double var;
    double ci_high;
    double ci_low;
    bool   informative;
    loc_t  type;

    dPrime() {
        this->dprime  = 0.0;
        this->chisq_p = false;
        this->var     = 0.0;
        this->ci_high = 0.0;
        this->ci_low  = 0.0;
        this->type    = uninformative;
    }
};

class PhasedSummary {
    map<string, int> sample_map;

public:
    uint     size;
    uint     sample_cnt;
    NucSum  *nucs;
    PhasedSample  *samples;
    dPrime **dprime;
    bool   **recomb;

    PhasedSummary(uint num_samples, uint num_genotypes) {
        this->sample_cnt = num_samples;
        this->samples    = new PhasedSample[this->sample_cnt];
        this->size       = num_genotypes;
        this->nucs       = new NucSum[this->size];
        this->dprime     = new dPrime *[this->size];
        for (uint i = 0; i < this->size; i++)
            this->dprime[i] = new dPrime[this->size];

        this->recomb     = new bool *[this->size];
        for (uint i = 0; i < this->size; i++) {
            this->recomb[i] = new bool[this->size];
            memset(this->recomb[i], 0, this->size);
        }
    }
    ~PhasedSummary() {
        if (this->nucs != NULL)
            delete [] this->nucs;
        if (this->dprime != NULL) {
            for (uint i = 0; i < this->size; i++)
                delete [] this->dprime[i];
            delete [] this->dprime;
        }
        if (this->recomb != NULL) {
            for (uint i = 0; i < this->size; i++)
                delete [] this->recomb[i];
            delete [] this->recomb;
        }
        if (this->samples != NULL)
            delete [] this->samples;
    }
    int add_sample(string name) {
        uint i = this->sample_map.size();
        this->sample_map[name] = i;
        this->samples[i].name = name;
        return i;
    }
};

class HBlock {
public:
    vector<int> loci;
    HBlock *next;

    HBlock() {
        this->next = NULL;
    }
};

class dPrimeBlocks {
    HBlock *_head;

public:
    dPrimeBlocks() {
        this->_head = NULL;
    }
    ~dPrimeBlocks() {
        HBlock *cur, *next;
        cur  = this->_head;
        next = cur->next;

        while (next != NULL) {
            delete cur;
            cur  = next;
            next = cur->next;
        }
    }
    HBlock *head() {
        return this->_head;
    }
    HBlock *initialize(set<int> &);
    HBlock *merge_adjacent(HBlock *);
    int     print();
};

void  help( void );
void  version( void );
int   parse_command_line(int, char**);
int   build_file_list(vector<pair<int, string> > &);
int   parse_population_map(string, map<string, int> &, map<int, int> &);

PhasedSummary *parse_fastphase(string);
PhasedSummary *parse_beagle(map<int, CSLocus *> &, string);
PhasedSummary *parse_beagle_haplotypes(map<int, CSLocus *> &, string);

int   summarize_phased_genotypes(PhasedSummary *);
int   calc_dprime(PhasedSummary *);
int   assign_alleles(NucSum, char &, char &, double &, double &);
int   write_dprime(string, PhasedSummary *);
int   four_gamete_test(string, map<string, int> &, PhasedSummary *, map<int, int> &, map<int, int> &);
int   dprime_blocks(string, map<string, int> &, PhasedSummary *, map<int, int> &, map<int, int> &);
bool  check_adjacent_blocks(PhasedSummary *, HBlock *);
int   enumerate_haplotypes(ofstream &, map<string, int> &, PhasedSummary *, uint, uint);
int   bucket_dprime(vector<double> &, vector<double> &, PhasedSummary *);
int   write_buckets(string, vector<double> &, vector<double> &);

#endif // __PHASEDSTACKS_H__
