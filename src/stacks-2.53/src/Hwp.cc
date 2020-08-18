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

#include "Hwp.h"

double
GuoThompson_Hwp::exec_locus(int start, int end, const Datum **d, size_t hap_cnt)
{
    std::ostream_iterator<unsigned long> out (cerr, " ");

    //
    // 1. Create several NxN diagonal matrix containing all N possible genotype combinations found in the population.
    //      g: the current matrix on the Markov chain.
    //      g_dswitch: the next potential matrix on the Markov chain if we d-switch.
    //      g_rswitch: the next potential matrix on the Markov chain if we r-switch.
    //
    HWMatrix g_1(hap_cnt), g_2(hap_cnt), g_3(hap_cnt);

    //
    // 2. Populate the matrix.
    //
    g_1.populate(start, end, d);

    // g_1.dump_matrix();
    // cerr << "Populated matrix. Number of alleles: " << g_1.alleles() << "; number of genotypes: " << g_1.genotypes() << "; number of hets: " << g_1.hets() << "\n";

    HWMatrix *g          = &g_1;
    HWMatrix *g_dswitch  = &g_2;
    HWMatrix *g_rswitch  = &g_3;
    double    log_rho    = 0.0;

    //
    // Markov chain burn-in period.
    //
    for (size_t n = 0; n < this->_burnin; n++)
        log_rho = walk_chain(&g, &g_dswitch, &g_rswitch, log_rho);
    
    double p_mean   = 0.0;
    double p_square = 0.0;
    double p_sim, se;

    //
    // Execute the Markov chain in b batches, with n steps in each batch. Average the estimated p-value.
    //
    for (size_t b = 0; b < this->_batches; b++) {

        size_t k = 0;

        for (size_t n = 0; n < this->_steps; n++) {
            log_rho = walk_chain(&g, &g_dswitch, &g_rswitch, log_rho);

            if (log_rho <= 0.0)
                k++;
        }

        p_sim     = (double) k / (double) this->_steps;
        p_mean   += p_sim;
        p_square += p_sim * p_sim;

        // cerr << "N: " << this->_steps << ", K: " << k << ", p-value: " << p_sim << "\n";
    }

    p_mean = p_mean / this->_batches;
    se     = p_square / ((double) this->_batches) / ((double) this->_batches - 1.0) - p_mean / ( (double) this->_batches - 1.0 ) * p_mean;
    se     = sqrt ( se );

    // cerr << "P-value: " << p_mean << " (" << se << " SE) \n";

    this->_p_value = p_mean;
    this->_se = se;

    return p_mean;
}

double
GuoThompson_Hwp::walk_chain(HWMatrix **g, HWMatrix **g_dswitch, HWMatrix **g_rswitch, double log_rho)
{
    HWMatrix *gptr;
    double    log_rho_prime = log_rho;

    //
    // 1. Randomly select two pairs of indexes to access two cells in the matrix.
    //    Make sure i_1 < i_2 and j_1 < j_2 to keep all valules in the lower diagonal of the matrix.
    //
    size_t i_1 = this->_indexes(this->_eng);
    size_t j_1 = this->_indexes(this->_eng);
    size_t i_2, j_2;
    do {
        i_2 = this->_indexes(this->_eng);
    } while (i_2 == i_1);
    do {
        j_2 = this->_indexes(this->_eng);
    } while (j_2 == j_1);

    if (i_2 < i_1)
        std::swap(i_1, i_2);
    if (j_2 < j_1)
        std::swap(j_1, j_2);
    
    //
    // 2. Calculate the value of delta and determine the switchability of the current matrix.
    //
    Switchable sw(*g, i_1, j_1, i_2, j_2);
    double pr_d, pr_r;

    if (sw.fully_switchable) {
        //
        // Propose both a new d-switched and r-switched matrix. 
        //
        pr_d = this->d_switch(**g, **g_dswitch, i_1, j_1, i_2, j_2, sw);
        pr_r = this->r_switch(**g, **g_rswitch, i_1, j_1, i_2, j_2, sw);

        //
        // Calculate the transition probability of the two switches based on equation 3 from Guo and Thompson, 1992.
        // We draw a random number between 0 and 1.0.
        //   If that number falls between 0 and pr_lim_d, we accept the dswitch transition.
        //   If that number falls between pr_lim_d and pr_lim_r, we accept the rswitch transition.
        //   Otherwise, we do not switch.
        //
        double pr_lim_d = this->transition_prob(pr_d);
        double pr_lim_r = pr_lim_d + this->transition_prob(pr_r);
        double tr_pr    = this->_transition(this->_eng);
            
        if (tr_pr <= pr_lim_d) {
            gptr           = *g;
            *g             = *g_dswitch;
            *g_dswitch     = gptr;
            log_rho_prime  = log_rho + log(pr_d);

        } else if (tr_pr <= pr_lim_r) {
            gptr           = *g;
            *g             = *g_rswitch;
            *g_rswitch     = gptr;
            log_rho_prime  = log_rho + log(pr_r);
        }

    } else if (sw.partially_switchable && sw.d_switchable) {
        //
        // Propose a new d-switched matrix and calculate the transistion probability.
        //
        pr_d = this->d_switch(**g, **g_dswitch, i_1, j_1, i_2, j_2, sw);

        //
        // Calculate the transition probability of switching based on equation 3 from Guo and Thompson, 1992.
        // We draw a random number between 0 and 1.0.
        //   If that number falls between 0 and pr_lim_d, we accept the dswitch transition.
        //   Otherwise, we do not switch.
        //
        double pr_lim_d = this->transition_prob(pr_d);
        double tr_pr    = this->_transition(this->_eng);

        if (tr_pr <= pr_lim_d) {
            gptr           = *g;
            *g             = *g_dswitch;
            *g_dswitch     = gptr;
            log_rho_prime  = log_rho + log(pr_d);
        }

    } else if (sw.partially_switchable && sw.r_switchable) {
        //
        // Propose a new r-switched matrix and calculate the transition probability.
        //
        pr_r = this->r_switch(**g, **g_rswitch, i_1, j_1, i_2, j_2, sw);

        //
        // Calculate the transition probability of the two switches based on equation 3 from Guo and Thompson, 1992.
        // We draw a random number between 0 and 1.0.
        //   If that number falls between 0 and pr_lim_r, we accept the rswitch transition.
        //   Otherwise, we do not switch.
        //
        double pr_lim_r = this->transition_prob(pr_r);
        double tr_pr    = this->_transition(this->_eng);

        if (tr_pr <= pr_lim_r) {
            gptr           = *g;
            *g             = *g_rswitch;
            *g_rswitch     = gptr;
            log_rho_prime  = log_rho + log(pr_r);
        }
    }
        
    return log_rho_prime;
}

inline double
GuoThompson_Hwp::d_switch(HWMatrix &g, HWMatrix &g_prime, size_t i_1, size_t j_1, size_t i_2, size_t j_2, Switchable &sw)
{
    //
    // Implementation of donation (d) switches, from Guo and Thompson, Table 1.
    //
    g_prime = g;

    if (sw.d0 || sw.d1) {
        g_prime.dec_genotype(i_1, j_1);
        g_prime.dec_genotype(i_2, j_2);
        g_prime.inc_genotype(i_1, j_2);
        g_prime.inc_genotype(i_2, j_1);
    } else if (sw.d2) {
        g_prime.dec_genotype(i_1, j_1);
        g_prime.dec_genotype(i_2, j_2);
        g_prime.inc_genotype(i_1, j_2);
        g_prime.inc_genotype(i_1, j_2);
    }

    double pr    = 0.0;
    double gamma = this->gamma(sw.Delta, i_1, j_1, i_2, j_2);

    if (sw.d0 || sw.d1) {
        pr = ( gamma * (double) g.cellv(i_1, j_1) * (double) g.cellv(i_2, j_2) ) / ( ((double) g.cellv(i_1, j_2) + 1.0) * ((double) g.cellv(i_2, j_1) + 1.0));

    } else if (sw.d2) {
        pr = ( gamma * (double) g.cellv(i_1, j_1) * (double) g.cellv(i_2, j_2) ) / ( ((double) g.cellv(i_1, j_2) + 2.0) * ((double) g.cellv(i_1, j_2) + 1.0));
    }

    return pr;
}


inline double
GuoThompson_Hwp::r_switch(HWMatrix &g, HWMatrix &g_prime, size_t i_1, size_t j_1, size_t i_2, size_t j_2, Switchable &sw)
{
    //
    // Implementation of reception (r) switches, from Guo and Thompson, Table 1.
    //
    g_prime = g;

    if (sw.r0 || sw.r1) {
        g_prime.inc_genotype(i_1, j_1);
        g_prime.inc_genotype(i_2, j_2);
        g_prime.dec_genotype(i_1, j_2);
        g_prime.dec_genotype(i_2, j_1);
    } else if (sw.r2) {
        g_prime.inc_genotype(i_1, j_1);
        g_prime.inc_genotype(i_2, j_2);
        g_prime.dec_genotype(i_1, j_2);
        g_prime.dec_genotype(i_1, j_2);
    }

    double pr    = 0.0;
    double gamma = this->gamma(sw.Delta, i_1, j_1, i_2, j_2);

    if (sw.r0 || sw.r1) {
        pr = ( (double) g.cellv(i_1, j_2) * (double) g.cellv(i_2, j_1)) / (gamma * ((double) g.cellv(i_1, j_1) + 1.0) * ((double) g.cellv(i_2, j_2) + 1.0));

    } else if (sw.r2) {
        pr = ( (double) g.cellv(i_1, j_2) * ((double) g.cellv(i_1, j_2) - 1.0) ) / (gamma * ((double) g.cellv(i_1, j_1) + 1.0) * ((double) g.cellv(i_2, j_2) + 1.0));
    }
    
    return pr;
}

inline double
GuoThompson_Hwp::transition_prob(double pr)
{
    //
    // Equation 3, Guo and Thompson 1992.
    //
    return std::min(1.0, pr) / 2.0;
}

inline double
GuoThompson_Hwp::gamma(size_t Delta, size_t i_1, size_t j_1, size_t i_2, size_t j_2)
{
    return (pow(2.0, Delta) * delta(i_1, j_1, i_2, j_2)) + (pow(2.0, (-1.0 * Delta)) * (1.0 - delta(i_1, j_1, i_2, j_2)));
}

inline double
GuoThompson_Hwp::delta(size_t i_1, size_t j_1, size_t i_2, size_t j_2)
{
    return (kronecker(i_1, j_1) + kronecker(i_2, j_2)) - (kronecker(i_1, j_1) * kronecker(i_2, j_2));
}

inline double
GuoThompson_Hwp::log_hwe_probability(HWMatrix &g)
{
    //
    // Implementation of equation 1 from Guo and Thompson, 1992.
    //
    size_t n_genotypes = 0;
    size_t hets        = 0;

    //
    // Sum up the total number of genotypes and the number of heterozygous individuals.
    //
    for (uint i = 0; i < g.alleles(); i++)
        for(uint j = 0; j <= i; j++) {
            if (j < i) hets += g.cellv(i, j);
            n_genotypes += g.cellv(i, j);
        }

    cerr << "Number of hets: " << hets << ", total number of genotypes: " << n_genotypes << "\n";
    cerr << "Number of hets: " << g.hets() << ", total number of genotypes: " << g.genotypes() << "\n";
    
    //
    // Calculate the nunmerator and denominator of the probability.
    //
    double num = 0.0;
    map<string, size_t> &hap_cnts = g.hap_cnts();
    for (map<string, size_t>::iterator i = hap_cnts.begin(); i != hap_cnts.end(); i++) {
        num += log_factorial(i->second);
    }

    double den = 0.0;
    for (uint i = 0; i < g.alleles(); i++)
        for(uint j = 0; j <= i; j++) {
            den += log_factorial(g.cellv(i, j));
        }

    double log_pr = ( (log_factorial(n_genotypes) + num) - (log_factorial(2*n_genotypes) + den) ) + (hets * log(2));

    return log_pr;
}

inline double
GuoThompson_Hwp::hwe_probability(HWMatrix &g)
{
    size_t n_genotypes = 0;
    size_t hets        = 0;

    //
    // Sum up the total number of genotypes and the number of heterozygous individuals.
    //
    for (uint i = 0; i < g.alleles(); i++)
        for(uint j = 0; j <= i; j++) {
            if (j < i) hets += g.cellv(i, j);
            n_genotypes += g.cellv(i, j);
        }

    double num = 1.0;
    map<string, size_t> &hap_cnts = g.hap_cnts();
    for (map<string, size_t>::iterator i = hap_cnts.begin(); i != hap_cnts.end(); i++) {
        num *= factorial(i->second);
    }
    
    double den = 1.0;
    for (uint i = 0; i < g.alleles(); i++)
        for(uint j = 0; j <= i; j++) {
            den *= factorial(g.cellv(i, j));
        }

    double pr = ( (factorial(n_genotypes) * num) / (factorial(2 * n_genotypes) * den) ) * pow(2, hets);

    return pr;
}

void
HWMatrix::dump_matrix()
{
    int  padding = 3;
    char buf[id_len];

    cerr << "Allele cnts:\n";
    for (map<string, size_t>::iterator it = this->_hap_cnts.begin(); it != this->_hap_cnts.end(); it++)
        cerr << it->first << ": " << it->second << "\n";
    
    cerr << "    ";
    for (uint i = 0; i < this->_hap_list.size(); i++)
        cerr << this->_hap_list[i] << "  ";
    cerr << "\n";

    for (uint i = 0; i < this->_hap_list.size(); i++) {
        cerr << this->_hap_list[i] << "  ";
    
        for(uint j = 0; j < this->_hap_list.size(); j++) {
            if (j <= i)
                snprintf(buf, id_len, "%*lu  ", padding, this->_g[i][j]);
            else
                snprintf(buf, id_len, "%*d  ", padding, 0);
            cerr << buf;
        }
        cerr << "\n";
    }
}

void
HWMatrix::populate(int start, int end, const Datum **d)
{
    for (int i = start; i <= end; i++) {
        if (d[i] == NULL)
            // No data, ignore this sample.
            continue;

        const vector<char*>& haps = d[i]->obshap;

        if (haps.size() > 2) {
            // Too many haplotypes, ignore this sample.
            continue;

        } else if (haps.size() == 1) {
            // Homozygote.
            if(!uncalled_haplotype(d[i]->obshap[0]))
                this->add_genotype(d[i]->obshap[0], d[i]->obshap[0]);

        } else {
            // Heterozygote.
            if(!uncalled_haplotype(d[i]->obshap[0]) &&
               !uncalled_haplotype(d[i]->obshap[1]))
                this->add_genotype(d[i]->obshap[0], d[i]->obshap[1]);
        }
    }
}

inline size_t
kronecker(size_t i, size_t j)
{
    return i == j ? 1 : 0;
}
