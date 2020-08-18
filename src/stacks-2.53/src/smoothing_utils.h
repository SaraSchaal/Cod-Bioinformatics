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

#ifndef __SMOOTHING_UTILS_H__
#define __SMOOTHING_UTILS_H__

#include <cmath>
#include <vector>

extern double sigma;

inline
double *
calc_weights()
{
    //
    // Calculate weights for window smoothing operations.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    //
    int     limit   = 3 * sigma;
    double *weights = new double[limit + 1];

    for (int i = 0; i <= limit; i++)
        weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(sigma, 2)));

    return weights;
}

template<class StatT>
int
determine_window_limits(vector<StatT *> &sites, uint center_bp, uint &pos_l, uint &pos_u)
{
    int limit   = 3 * sigma;
    int limit_l = center_bp - limit > 0 ? center_bp - limit : 0;
    int limit_u = center_bp + limit;

    while (pos_l < sites.size()) {
        if (sites[pos_l] == NULL) {
            pos_l++;
        } else {
            if (sites[pos_l]->bp < limit_l)
                pos_l++;
            else
                break;
        }
    }
    while (pos_u < sites.size()) {
        if (sites[pos_u] == NULL) {
            pos_u++;
        } else {
            if (sites[pos_u]->bp < limit_u)
                pos_u++;
            else
                break;
        }
    }
    return 0;
}

#endif // __SMOOTHING_UTILS_H__
