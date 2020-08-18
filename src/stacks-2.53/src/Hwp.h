// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2018, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __HWP_H__
#define __HWP_H__

#include <random>
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::random_device;
using std::mt19937;
using std::seed_seq;

#include <iterator>
#include "constants.h"
#include "utils.h"
#include "PopMap.h"
#include "PopSum.h"

//
// Utility identity function used in both GuoThompson_Hwp and Switchable.
//
size_t kronecker(size_t i, size_t j);

class HWMatrix {
    vector<string>      _hap_list;
    map<string, size_t> _hap_index;
    map<size_t, string> _rev_index;
    map<string, size_t> _hap_cnts;
    size_t   _hets;
    size_t   _n_genotypes;
    size_t   _n_alleles;
    size_t **_g;

 public:
    HWMatrix(size_t n_alleles) {
        this->_n_alleles   = n_alleles;
        this->_n_genotypes = 0;
        this->_hets        = 0;
        this->_g           = new size_t * [this->_n_alleles];
        for (uint i = 0; i < this->_n_alleles; i++) {
            this->_g[i] = new size_t[this->_n_alleles];
            memset(this->_g[i], 0, this->_n_alleles * sizeof(size_t));
        }
    }
    HWMatrix(const HWMatrix &rhs) {
        this->_n_alleles   = rhs._n_alleles;
        this->_n_genotypes = rhs._n_genotypes;
        this->_hets        = rhs._hets;
        this->_g           = new size_t * [this->_n_alleles];

        for (uint i = 0; i < this->_n_alleles; i++) {
            this->_g[i] = new size_t[this->_n_alleles];
            memset(this->_g[i], 0, this->_n_alleles * sizeof(size_t));
            for (uint j = 0; j <= i; j++)
                this->_g[i][j] = rhs._g[i][j];
        }
        this->_hap_list  = rhs._hap_list;
        this->_hap_index = rhs._hap_index;
        this->_rev_index = rhs._rev_index;
        this->_hap_cnts  = rhs._hap_cnts;
    }
    ~HWMatrix() {
        for (uint i = 0; i < this->_n_alleles; i++)
            delete [] this->_g[i];
        delete [] this->_g;
    }
    HWMatrix &operator= (const HWMatrix &rhs) {
        assert(this->_n_alleles == rhs._n_alleles);

        this->_n_alleles   = rhs._n_alleles;
        this->_n_genotypes = rhs._n_genotypes;
        this->_hets        = rhs._hets;

        for (uint i = 0; i < this->_n_alleles; i++)
            for (uint j = 0; j <= i; j++) {
                this->_g[i][j] = rhs._g[i][j];
            }

        this->_hap_list  = rhs._hap_list;
        this->_hap_index = rhs._hap_index;
        this->_rev_index = rhs._rev_index;
        this->_hap_cnts  = rhs._hap_cnts;

        return *this;
    }
    size_t &cellr(size_t i, size_t j) {
        //
        // Only process values in the lower triangular matrix.
        //
        return (j > i) ? this->_g[j][i] : this->_g[i][j];
    }
    size_t cellv(size_t i, size_t j) const {
        //
        // Only process values in the lower triangular matrix.
        //
        return (j > i) ? this->_g[j][i] : this->_g[i][j];
    }
    size_t alleles() const { return this->_n_alleles; }
    size_t hets() const { return this->_hets; }
    size_t genotypes() const { return this->_n_genotypes; }
    map<string, size_t>& hap_cnts() { return this->_hap_cnts; }

    void add_genotype(string allele_1, string allele_2) {
        if (this->_hap_index.count(allele_1) == 0) {
            _hap_index[allele_1] = _hap_list.size();
            _rev_index[_hap_list.size()] = allele_1;
            _hap_cnts[allele_1]  = 1;
            _hap_list.push_back(allele_1);        
        } else {
            _hap_cnts[allele_1]++;
        }
        if (this->_hap_index.count(allele_2) == 0) {
            _hap_index[allele_2] = _hap_list.size();
            _rev_index[_hap_list.size()] = allele_2;
            _hap_cnts[allele_2]  = 1;
            _hap_list.push_back(allele_2);        
        } else {
            _hap_cnts[allele_2]++;
        }
        
        size_t i = this->_hap_index[allele_1];
        size_t j = this->_hap_index[allele_2];
        
        if (j > i)
            this->_g[j][i]++;
        else
            this->_g[i][j]++;

        if (i != j)
            this->_hets++;

        this->_n_genotypes++;
    }
    size_t inc_genotype(size_t i, size_t j) {
        assert(i < this->_n_alleles && j < this->_n_alleles);
        
        if (j > i)
            this->_g[j][i]++;
        else
            this->_g[i][j]++;

        if (i != j) this->_hets++;
        
        this->_n_genotypes++;

        return (j > i) ? this->_g[j][i] : this->_g[i][j];
    }
    size_t dec_genotype(size_t i, size_t j) {
        assert(i < this->_n_alleles && j < this->_n_alleles);

        if (j > i)
            this->_g[j][i]--;
        else
            this->_g[i][j]--;

        if (i != j) this->_hets--;

        this->_n_genotypes--;

        return (j > i) ? this->_g[j][i] : this->_g[i][j];
    }
    void populate(int, int, const Datum **);
    void dump_matrix();
};

class Switchable {
public:
    size_t Delta;
    bool d0;    // Type 0 donation switch, AKA d0-switchable.
    bool r0;    // Type 0 reception switch, AKA r0-switchable.
    bool d1;    // Type 1 donation switch, AKA d1-switchable.
    bool r1;    // Type 1 reception switch, AKA r1-switchable.
    bool d2;    // Type 2 donation switch, AKA d2-switchable.
    bool r2;    // Type 2 reception switch, AKA r2-switchable.
    bool d_switchable;
    bool r_switchable;
    bool partially_switchable;
    bool fully_switchable;

    Switchable() {
        this->Delta = 0;
        this->d0 = false;
        this->r0 = false;
        this->d1 = false;
        this->r1 = false;
        this->d2 = false;
        this->r2 = false;
        this->d_switchable = false;
        this->r_switchable = false;
        this->partially_switchable = false;
        this->fully_switchable     = false;
    }
    Switchable(HWMatrix *g, size_t i_1, size_t j_1, size_t i_2, size_t j_2) {
        this->Delta = 0;
        this->d0 = false;
        this->r0 = false;
        this->d1 = false;
        this->r1 = false;
        this->d2 = false;
        this->r2 = false;
        this->d_switchable = false;
        this->r_switchable = false;
        this->partially_switchable = false;
        this->fully_switchable     = false;

        this->switchability(g, i_1, j_1, i_2, j_2);
    }

    size_t switchability(HWMatrix *g, size_t i_1, size_t j_1, size_t i_2, size_t j_2) {
        this->Delta = kronecker(i_1, j_1) + kronecker(i_1, j_2) + kronecker(i_2, j_1) + kronecker(i_2, j_2);

        assert(this->Delta >= 0 && this->Delta <= 2);

        //
        // Determine if and how this set of indices are switchable.
        //
        // Delta can be 0, 1, or 2 only. If the right cells in the contingency table have counts > 0
        // then the state is switchabe, either a donation or reception switch (paritally switchable)
        // or both (fully switchable).
        //
        switch(Delta) {
        case 0:
            if (g->cellv(i_1, j_1) * g->cellv(i_2, j_2) != 0) d0 = true;
            if (g->cellv(i_1, j_2) * g->cellv(i_2, j_1) != 0) r0 = true;
            break;
        case 1:
            if (g->cellv(i_1, j_1) * g->cellv(i_2, j_2) != 0) d1 = true;
            if (g->cellv(i_1, j_2) * g->cellv(i_2, j_1) != 0) r1 = true; 
            break;
        case 2:
            if (g->cellv(i_1, j_1) * g->cellv(i_2, j_2) != 0) d2 = true;
            if (g->cellv(i_1, j_2) >= 2)                      r2 = true;
            break;
        }

        this->d_switchable = (this->d0 || this->d1 || this->d2) ? true : false;
        this->r_switchable = (this->r0 || this->r1 || this->r2) ? true : false;
        
        if (this->d_switchable && this->r_switchable) {
            this->fully_switchable = true;

        } else if (this->d_switchable || this->r_switchable) {
            this->partially_switchable = true;
        }
        // else: not switchable.

        return this->Delta;
    }
};

class GuoThompson_Hwp {
    const size_t _burnin  = 10000;
    const size_t _steps   = 10000;
    const size_t _batches = 20;

    random_device _r;
    mt19937       _eng;
    //
    // Random number generator for choosing whether to transition between states.
    //
    uniform_real_distribution<double> _transition;
    //
    // Random number generator for choosing cells from our genotype matrix.
    //
    uniform_int_distribution<uint16_t> _indexes;

public:
    size_t _n_alleles;
    double _p_value;
    double _se;

    GuoThompson_Hwp(size_t n_alleles) {
        this->_n_alleles  = n_alleles;

        seed_seq seed{_r(), _r(), _r(), _r(), _r(), _r(), _r(), _r()};
        _eng.seed(seed);

        this->_indexes    = uniform_int_distribution<uint16_t>(0, n_alleles - 1);
        this->_transition = uniform_real_distribution<double>(0, 1);
    }

    double exec_locus(int, int, const Datum **, size_t);
    double log_hwe_probability(HWMatrix &);
    double hwe_probability(HWMatrix &);

private:
    double walk_chain(HWMatrix **, HWMatrix **, HWMatrix **, double);
    double r_switch(HWMatrix &, HWMatrix &, size_t i_1, size_t j_1, size_t i_2, size_t j_2, Switchable &sw);
    double d_switch(HWMatrix &, HWMatrix &, size_t i_1, size_t j_1, size_t i_2, size_t j_2, Switchable &sw);
    double delta(size_t i_1, size_t j_1, size_t i_2, size_t j_2);
    double gamma(size_t Delta, size_t i_1, size_t j_1, size_t i_2, size_t j_2);
    double transition_prob(double pr);
};

#endif // __HWP_H__
