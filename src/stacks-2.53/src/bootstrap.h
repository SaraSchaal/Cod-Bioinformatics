// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2014, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__

#include <cmath>
#include <vector>

#include "smoothing_utils.h"

extern double   sigma;
extern int      bootstrap_reps;
extern bool     bootstrap_wl;
extern set<int> bootstraplist;

//
// Bootstrap resamplign structure.
//
class BSample {
public:
    int    bp;
    int    alleles;
    bool   fixed;
    double stat[PopStatSize];

    BSample() {
        this->bp      = 0;
        this->alleles = 0;
        this->fixed   = false;
        for (uint i = 0; i < PopStatSize; i++)
            this->stat[i] = 0.0;

    }
};

template<class StatT=PopStat>
class Bootstrap {
    double                 *weights; // Weight matrix to apply while smoothing.
    vector<vector<double> > stats;
    uint                    num_stats;

public:
    Bootstrap(uint size)  {
        this->num_stats = size;
        this->weights   = calc_weights();
        this->stats.resize(size, vector<double>());
    }
    ~Bootstrap() {
        delete [] this->weights;
    }

    int    add_data(vector<StatT *> &);
    int    execute(vector<StatT *> &);
    int    execute_mixed(vector<StatT *> &);
    double pval(double, vector<double> &);
};

template<class StatT>
int
Bootstrap<StatT>::add_data(vector<StatT *> &sites)
{
    for (uint i = 0; i < sites.size(); i++) {
        if (sites[i] != NULL && sites[i]->fixed == false)
            for (uint j = 0; j < this->num_stats; j++)
                this->stats[j].push_back(sites[i]->stat[j]);
    }

    return 0;
}

template<class StatT>
int
Bootstrap<StatT>::execute(vector<StatT *> &sites)
{
    #pragma omp parallel
    {
        PopStat *c;
        double final_weight, sum, weighted_stat[PopStatSize];
        int  dist, index;
        uint pos_l = 0;
        uint pos_u = 0;

        //#pragma omp for schedule(dynamic, 1)
        for (uint pos_c = 0; pos_c < sites.size(); pos_c++) {
            c = sites[pos_c];

            if (c == NULL)
                continue;

            if (bootstrap_wl && bootstraplist.count(c->loc_id) == 0)
                continue;

            // cerr << "Bootstrapping " << c->loc_id << "; pos_c: " << pos_c << "; bp: " << c->bp << "\n";

            determine_window_limits(sites, c->bp, pos_l, pos_u);

            int size = 0;
            for (uint i = pos_l; i < pos_u;  i++)
                if (sites[i] != NULL) size++;

            //
            // Allocate an array of bootstrap resampling objects.
            //
            BSample *bs = new BSample[size];

            //
            // Populate the BSample objects.
            //
            int j = 0;
            for (uint i = pos_l; i < pos_u;  i++) {
                if (sites[i] == NULL) continue;
                bs[j].bp      = sites[i]->bp;
                bs[j].alleles = sites[i]->alleles;
                j++;
            }

            vector<vector<double> > resampled_stats(this->num_stats, vector<double>());
            for (uint i = 0; i < this->num_stats; i++)
                resampled_stats[i].reserve(bootstrap_reps);

            //
            // Bootstrap this bitch.
            //
            for (int i = 0; i < bootstrap_reps; i++) {
                // if (i % 100 == 0) cerr << "      Bootsrap rep " << i << "\n";

                for (uint k = 0; k < this->num_stats; k++)
                    weighted_stat[k] = 0.0;
                sum = 0.0;

                for (j = 0; j < size; j++) {
                    //
                    // Distance from center of window.
                    //
                    dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;
                    //
                    // Resample for this round of bootstrapping.
                    //
                    index = (int) (this->stats[0].size() * (random() / (RAND_MAX + 1.0)));
                    for (uint k = 0; k < this->num_stats; k++)
                        bs[j].stat[k] = this->stats[k][index];

                    final_weight = (bs[j].alleles - 1) * this->weights[dist];
                    for (uint k = 0; k < this->num_stats; k++)
                        weighted_stat[k] += bs[j].stat[k] * final_weight;
                    sum += final_weight;
                }

                // cerr << "    New weighted Fst value: " << weighted_fst / sum << "\n";
                for (uint k = 0; k < this->num_stats; k++)
                    resampled_stats[k].push_back(weighted_stat[k] / sum);
            }

            //
            // Cacluate the p-value for this window based on the empirical Fst distribution.
            //
            for (uint k = 0; k < this->num_stats; k++) {
                sort(resampled_stats[k].begin(), resampled_stats[k].end());
                c->bs[k] = this->pval(c->smoothed[k], resampled_stats[k]);
            }

            delete [] bs;
        }
    }

    return 0;
}

template<class StatT>
int
Bootstrap<StatT>::execute_mixed(vector<StatT *> &sites)
{
    #pragma omp parallel
    {
        PopStat *c;
        double final_weight, sum, weighted_stat[PopStatSize];
        int  dist, index;
        uint pos_l = 0;
        uint pos_u = 0;

        //#pragma omp for schedule(dynamic, 1)
        for (uint pos_c = 0; pos_c < sites.size(); pos_c++) {
            c = sites[pos_c];

            if (c == NULL || c->fixed == true)
                continue;

            if (bootstrap_wl && bootstraplist.count(c->loc_id) == 0)
                continue;

            // cerr << "Bootstrapping " << c->loc_id << "; pos_c: " << pos_c << "; bp: " << c->bp << "\n";

            determine_window_limits(sites, c->bp, pos_l, pos_u);

            int size = 0;
            for (uint i = pos_l; i < pos_u;  i++)
                if (sites[i] != NULL) size++;

            //
            // Allocate an array of bootstrap resampling objects.
            //
            BSample *bs = new BSample[size];

            //
            // Populate the BSample objects.
            //
            int j = 0;
            for (uint i = pos_l; i < pos_u;  i++) {
                if (sites[i] == NULL)
                    continue;
                bs[j].bp      = sites[i]->bp;
                bs[j].alleles = sites[i]->alleles;
                bs[j].fixed   = sites[i]->fixed;
                for (uint k = 0; k < this->num_stats; k++)
                    bs[j].stat[k] = sites[i]->stat[k];
                j++;
            }

            //
            // Precompute the fraction of the window that will not change during resampling.
            //
            double partial_weighted_stat[this->num_stats];
            double partial_sum = 0.0;
            memset(partial_weighted_stat, 0, this->num_stats);

            for (j = 0; j < size; j++) {
                if (bs[j].fixed == false) continue;

                dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

                final_weight  = (bs[j].alleles - 1.0) * this->weights[dist];
                partial_sum  += final_weight;
                for (uint k = 0; k < this->num_stats; k++)
                    partial_weighted_stat[k] += bs[j].stat[k] * final_weight;
            }

            vector<vector<double> > resampled_stats(this->num_stats, vector<double>());
            for (uint i = 0; i < this->num_stats; i++)
                resampled_stats[i].reserve(bootstrap_reps);

            // cerr << "Window starts at " << bs[0].bp << "; centered on " << c->bp << "\n";

            //
            // Bootstrap this bitch.
            //
            for (int i = 0; i < bootstrap_reps; i++) {
                // if (i % 100 == 0) cerr << "      Bootsrap rep " << i << "\n";

                for (uint k = 0; k < this->num_stats; k++)
                    weighted_stat[k] = partial_weighted_stat[k];
                sum = partial_sum;

                for (j = 0; j < size; j++) {
                    if (bs[j].fixed == true) continue;

                    dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

                    //
                    // Resample for this round of bootstrapping.
                    //
                    index = (int) (this->stats[0].size() * (random() / (RAND_MAX + 1.0)));
                    for (uint k = 0; k < this->num_stats; k++)
                        bs[j].stat[k] = this->stats[k][index];

                    final_weight = (bs[j].alleles - 1) * this->weights[dist];
                    for (uint k = 0; k < this->num_stats; k++)
                        weighted_stat[k] += bs[j].stat[k] * final_weight;
                    sum += final_weight;
                }

                // cerr << "    New weighted value: " << (weighted_stat[0] / sum) << "\n";
                for (uint k = 0; k < this->num_stats; k++)
                    resampled_stats[k].push_back(weighted_stat[k] / sum);
            }

            //
            // Cacluate the p-value for this window based on the empirical Fst distribution.
            //
            for (uint k = 0; k < this->num_stats; k++) {
                sort(resampled_stats[k].begin(), resampled_stats[k].end());
                c->bs[k] = this->pval(c->smoothed[k], resampled_stats[k]);
            }

            delete [] bs;
        }
    }

    return 0;
}

template<class StatT>
double
Bootstrap<StatT>::pval(double stat, vector<double> &dist)
{
    // Rank `stat` relative to the bootstrapped distribution.
    // `first_greater` is between 0 (if `stat` is smaller than all values in the
    // distribution) and `dist.size()` (if `stat` is greater than all of them).
    size_t first_greater = upper_bound(dist.begin(), dist.end(), stat) - dist.begin();

    // cerr << "Generated Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    //         cerr << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cerr << "Comparing Fst value: " << stat
    //          << " at position " << (up - dist.begin()) << " out of "
    //          << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return 1.0 - (double) first_greater / dist.size();
}

#endif // __BOOTSTRAP_H__
