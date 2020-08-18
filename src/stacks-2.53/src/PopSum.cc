// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2018, Julian Catchen <jcatchen@illinois.edu>
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

#include "PopSum.h"

LocPopSum::LocPopSum(size_t cloc_len, const MetaPopInfo& mpopi)
{
    this->_pop_cnt  = mpopi.pops().size();
    this->_meta_pop = new LocTally(cloc_len);
    this->_per_pop  = new LocSum * [this->_pop_cnt];
    this->_hapstats_per_pop = new LocStat * [this->_pop_cnt];

    for (size_t i = 0; i < mpopi.pops().size(); i++) {
        this->_per_pop[i] = new LocSum(cloc_len);
        this->_hapstats_per_pop[i] = NULL;
    }
}

LocPopSum::~LocPopSum()
{
    for (uint i = 0; i < this->_pop_cnt; i++) {
        delete this->_per_pop[i];
        delete this->_hapstats_per_pop[i];
    }
    delete [] this->_per_pop;
    delete [] this->_hapstats_per_pop;
    delete this->_meta_pop;
}

int
LocPopSum::sum_pops(const CSLocus *cloc, Datum const*const* d, const MetaPopInfo &mpopi,
                    bool verbose, ostream &log_fh)
{
    uint len = strlen(cloc->con);
    int res;
    set<int> snp_cols;
    int incompatible_loci = 0;

    for (uint i = 0; i < this->_pop_cnt; i++) {
        LocSum *s = this->_per_pop[i];

        const Pop& pop = mpopi.pops().at(i);
        //
        // Check if this locus has already been filtered and is NULL in all individuals.
        //
        bool filtered = true;
        for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
            if (d[k] != NULL)
                filtered = false;
        }
        if (filtered == true) {
            for (uint k = 0; k < len; k++) {
                s->nucs[k].filtered_site = true;
            }
            continue;
        }

        //
        // The catalog records which nucleotides are heterozygous. For these nucleotides we will
        // calculate observed genotype frequencies, allele frequencies, and expected genotype frequencies.
        //
        for (uint k = 0; k < cloc->snps.size(); k++) {
            res = this->tally_heterozygous_pos(cloc, d, s,
                                               cloc->snps[k]->col, k, pop.first_sample, pop.last_sample);
            //
            // If site is incompatible (too many alleles present), log it.
            //
            if (res < 0) {
                s->nucs[cloc->snps[k]->col].incompatible_site = true;
                incompatible_loci++;
                if (verbose)
                    log_fh << "within_population\t"
                           << "incompatible_locus\t"
                           << cloc->id << "\t"
                           << cloc->loc.chr() << "\t"
                           << cloc->sort_bp(cloc->snps[k]->col) +1 << "\t"
                           << cloc->snps[k]->col << "\t"
                           << pop.name << "\n";
                DOES_NOT_HAPPEN;
                // @nick 2018-12-17: I don't think this should ever happen with genotypes
                // provided by gstacks, or generally, loaded from a VCF (as diploidy is
                // implicit in the format).
            }

            snp_cols.insert(cloc->snps[k]->col);
        }
        //
        // For all other fixed sites, we just need to record them.
        //
        for (uint k = 0; k < len; k++) {
            if (snp_cols.count(k)) continue;
            this->tally_fixed_pos(cloc, d, s, k, pop.first_sample, pop.last_sample);
        }

        snp_cols.clear();
    }

    return 0;
}

int
LocPopSum::tally_fixed_pos(const CSLocus *cloc, Datum const*const* d, LocSum *s,
                           int pos, uint start, uint end)
{
    double num_indv = 0.0;
    char   p_nuc    = 0;

    s->nucs[pos].reset();

    for (uint i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len)
            continue;
        //
        // Before counting this individual, make sure the model definitively called this
        // position as hEterozygous or hOmozygous.
        //
        assert(d[i]->model[pos] != 'E');
        num_indv++;
        p_nuc = cloc->con[pos];
    }
    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = cloc->id;
    s->nucs[pos].bp       = cloc->sort_bp(pos);
    s->nucs[pos].fixed    = true;
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = 2 * num_indv;

    if (num_indv > 0) {
        s->nucs[pos].p        =  1.0;
        s->nucs[pos].p_nuc    =  p_nuc;
        s->nucs[pos].obs_hom  =  1.0;
        s->nucs[pos].obs_het  =  0.0;
        s->nucs[pos].exp_hom  =  1.0;
        s->nucs[pos].exp_het  =  0.0;
        s->nucs[pos].stat[0]  =  0.0; // pi
        s->nucs[pos].stat[1]  = -7.0; // fis
        s->nucs[pos].stat[2]  =  0.0; // HWE deviation
    }

    return 0;
}

int
LocPopSum::tally_heterozygous_pos(const CSLocus *cloc, Datum const*const* d, LocSum *s,
                                  int pos, int snp_index, uint start, uint end)
{
    //
    // Tally up the genotype frequencies.
    //
    int  nucs[4] = {0};
    uint i;
    char nuc;

    s->nucs[pos].reset();
    //cout << "  Calculating summary stats at het locus " << cloc->id << " position " << pos << "; snp_index: " << snp_index << "\n";

    //
    // Iterate over each individual in this sub-population.
    //
    for (i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len || d[i]->model[pos] == 'U') continue;

        //
        // Pull each allele for this SNP from the observed haplotype.
        //
        for (uint j = 0; j < d[i]->obshap.size(); j++) {
            nuc = d[i]->obshap[j][snp_index];

            switch(nuc) {
            case 'A':
            case 'a':
                nucs[0]++;
                break;
            case 'C':
            case 'c':
                nucs[1]++;
                break;
            case 'G':
            case 'g':
                nucs[2]++;
                break;
            case 'T':
            case 't':
                nucs[3]++;
                break;
            }
        }
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2)
        return -1;

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    char p_allele = 0;
    char q_allele = 0;

    i = 0;
    while (p_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 0:
                p_allele = 'A';
                break;
            case 1:
                p_allele = 'C';
                break;
            case 2:
                p_allele = 'G';
                break;
            case 3:
                p_allele = 'T';
                break;
            }
        }
        i++;
    }
    while (q_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 1:
                q_allele = 'C';
                break;
            case 2:
                q_allele = 'G';
                break;
            case 3:
                q_allele = 'T';
                break;
            }
        }
        i++;
    }
    //cout << "  P Allele: " << p_allele << "; Q Allele: " << q_allele << "\n";

    //
    // Calculate observed genotype frequencies.
    //
    double num_indv = 0.0;
    double obs_het  = 0.0;
    double obs_p    = 0.0;
    double obs_q    = 0.0;

    for (i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len) continue;
        //
        // Before counting this individual, make sure the model definitively called this
        // position as hEterozygous or hOmozygous.
        //
        if (d[i]->model[pos] == 'E' || d[i]->model[pos] == 'O')
            num_indv++;
        else
            continue;

        if (d[i]->obshap.size() > 1 &&
            this->tally_observed_haplotypes(d[i]->obshap, snp_index) == 2)
                obs_het++;
        else if (d[i]->obshap[0][snp_index] == p_allele)
            obs_p++;
        else if (d[i]->obshap[0][snp_index] == q_allele)
            obs_q++;
    }
    //cout << "  Num Individuals: " << num_indv << "; Obs Hets: " << obs_het << "; Obs P: " << obs_p << "; Obs Q: " << obs_q << "\n";

    if (num_indv == 0) return 0;

    //
    // Calculate total number of alleles
    //
    double tot_alleles = num_indv * 2;
    double allele_p    = obs_het + (2 * obs_p);
    double allele_q    = obs_het + (2 * obs_q);

    //
    // Calculate Pi, equivalent to expected heterozygosity (exp_het)
    //
    s->nucs[pos].stat[0] = this->pi(tot_alleles, allele_p, allele_q);

    if (s->nucs[pos].stat[0] == 0.0)
        s->nucs[pos].fixed = true;

    double hwe_pval = 0.0;
    //
    // Calculates deviation from Hardy-Weinberg equilibrium.
    //
    if (calc_hwp) {
        hwe_pval = this->hwe(num_indv, allele_p, allele_q, obs_p, obs_q, obs_het);

        s->nucs[pos].stat[2] = hwe_pval;
    }

    //
    // Convert to allele frequencies
    //
    allele_p = allele_p / tot_alleles;
    allele_q = allele_q / tot_alleles;

    //cout << "  P allele frequency: " << allele_p << "; Q allele frequency: " << allele_q << "\n";

    //
    // Calculate expected genotype frequencies.
    //
    double exp_het = 2 * allele_p * allele_q; // 2pq
    // double exp_p   = allele_p * allele_p;     // p^2
    // double exp_q   = allele_q * allele_q;     // q^2

    //cout << "  Expected Het: " << exp_het << "; Expected P: " << exp_p << "; Expected Q: " << exp_q << "\n";

    obs_het = obs_het / num_indv;
    obs_p   = obs_p   / num_indv;
    obs_q   = obs_q   / num_indv;

    //cout << "  Obs Hets Freq: " << obs_het << "; Obs P Freq: " << obs_p << "; Obs Q Freq: " << obs_q << "\n";

    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = cloc->id;
    s->nucs[pos].bp       = cloc->sort_bp(pos);
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = tot_alleles;
    s->nucs[pos].p        = allele_p > allele_q ? allele_p : allele_q;
    s->nucs[pos].p_nuc    = allele_p > allele_q ? p_allele : q_allele;
    s->nucs[pos].q_nuc    = allele_p > allele_q ? q_allele : p_allele;
    s->nucs[pos].obs_hom  = 1 - obs_het;
    s->nucs[pos].obs_het  = obs_het;
    s->nucs[pos].exp_hom  = 1 - exp_het;
    s->nucs[pos].exp_het  = exp_het;

    //
    // Calculate F_is, the inbreeding coefficient of an individual (I) relative to the subpopulation (S).
    // F_is measures reductions in observed heterozygosity with respect to that expected under Hardy-Weinberg
    // Equilibrium.
    //   Fis = (exp_het - obs_het) / exp_het
    //
    double fis = s->nucs[pos].pi == 0 ? -7 : (s->nucs[pos].pi - obs_het) / s->nucs[pos].pi;

    s->nucs[pos].stat[1] = fis;

    return 0;
}

double
LocPopSum::pi(double tot_alleles, double p, double q)
{
    //
    // Calculate Pi, equivalent to expected heterozygosity:
    //  pi = 1 - Sum_i( (n_i choose 2) ) / (n choose 2)
    //
    double pi =
        LocPopSum::binomial_coeff(p, 2) +
        LocPopSum::binomial_coeff(q, 2);
    pi = pi / LocPopSum::binomial_coeff(tot_alleles, 2);
    pi = 1 - pi;

    return pi;
}

double
LocPopSum::binomial_coeff(double n, double k)
{
    if (n < k) return 0.0;
    //
    // Compute the binomial coefficient using the method of:
    // Y. Manolopoulos, "Binomial coefficient computation: recursion or iteration?",
    // ACM SIGCSE Bulletin, 34(4):65-67, 2002.
    //
    double r = 1.0;
    double s = (k < n - k) ? n - k + 1 : k + 1;

    for (double i = n; i >= s; i--)
        r = r * i / (n - i + 1);

    return r;
}

double
LocPopSum::hwe(double n, double p, double q, double p_hom, double q_hom, double hets)
{
    //
    // Compute the deviation from Hardy-Weinberg equilibrium using an exact test.
    //   Weir, Genetic Data Analysis II, 1996, Chapter 3, pp98-100.
    //   Wray and Visscher. Population genetics and its relevance to gene mapping. Chapter 6 in
    //     Statistical Genetics: Gene Mapping Through Linkage and Association. pp 90-91.
    //

    double p_fac = 0.0;
    for (uint i = 1; i <= p; i++) p_fac += log(i);
    double q_fac = 0.0;
    for (uint i = 1; i <= q; i++) q_fac += log(i);
    double n_fac = 0.0;
    for (uint i = n+1; i <= 2*n; i++) n_fac += log(i);

    //
    // To get the exact value we must also calculate all the probabilities that are less likely
    // than the observed probability. That is, given the same marginal totals, those with fewer heterozygotes.
    // To do this we will subtract 2 from hets and add 1 each to the homozygote classes until
    // we reach 1 or 0 hets.
    //
    // Defined in:
    //   Louis and Dempster, 1987. An Exact Test for Hardy-Weinberg and Multiple Alleles. Biometrics, Vol. 43, No. 4.
    //
    double s_hets, s_p_hom, s_q_hom;

    if (p <= q) {
        s_p_hom = 0;
        s_hets  = p;
        s_q_hom = (q - p) / 2.0;
    } else {
        s_p_hom = (p - q) / 2.0;
        s_hets  = q;
        s_q_hom = 0;
    }

    double obs_log_pr = log_hwp_pr(n_fac, p_fac, q_fac, p_hom, q_hom, hets);
    double hwe_pr     = 0.0;

    while (s_hets >= 0) {
        double log_pr = log_hwp_pr(n_fac, p_fac, q_fac, s_p_hom, s_q_hom, s_hets);
        double pr     = exp(log_pr);
        hwe_pr += log_pr > obs_log_pr ? 0.0 : pr;

        s_hets  -= 2;
        s_p_hom += 1;
        s_q_hom += 1;
    }

    return hwe_pr;
}

inline double
LocPopSum::log_hwp_pr(double n_fac, double p_fac, double q_fac, double pp, double qq, double pq)
{
    //
    // [n! (n_p)! (n_q)! 2^(n_pq)] / [(n_pp)! (n_pq)! (n_qq)! (2n)! ]
    //
    double num = p_fac + q_fac + (pq * log(2));

    double a = 0.0;
    for (uint i = 1; i <= pp; i++)  a += log(i);
    double b = 0.0;
    for (uint i = 1; i <= pq; i++)  b += log(i);
    double c = 0.0;
    for (uint i = 1; i <= qq; i++)  c += log(i);
    double den = a + b + c + n_fac;

    // cerr << " n_fac: " << n_fac
    //      << " p_fac: " << p_fac
    //      << " q_fac: " << q_fac
    //      << " pp: " << pp
    //      << " qq: " << qq
    //      << " pq: " << pq
    //      << " pp_fac: " << a
    //      << " pq_fac: " << b
    //      << " qq_fac: " << c
    //      << " num: " << num
    //      << " den: " << den
    //      << " log(hwe_pr): " << num - den
    //      << " hwe_pr: " << exp(num - den) << "\n";

    return num - den;
}

int
LocPopSum::tally_metapop(const CSLocus *cloc)
{
    int       variable_pop;
    uint16_t  p_cnt, q_cnt, col;
    LocSum  **s  = this->_per_pop;
    LocTally *mp = this->_meta_pop;

    for (col = 0; col < cloc->len; col++) {

        mp->nucs[col].reset();

        mp->nucs[col].col    = col;
        mp->nucs[col].bp     = cloc->sort_bp(col);
        mp->nucs[col].loc_id = cloc->id;

        this->tally_ref_alleles(col,
                                mp->nucs[col].allele_cnt,
                                mp->nucs[col].p_allele,
                                mp->nucs[col].q_allele,
                                p_cnt, q_cnt);

        //
        // Is this site variable?
        //
        if (mp->nucs[col].allele_cnt > 1)
            mp->nucs[col].fixed = false;

        for (uint j = 0; j < this->_pop_cnt; j++) {
            //
            // Sum the number of individuals examined at this locus across populations.
            //
            mp->nucs[col].num_indv += s[j]->nucs[col].num_indv;
            mp->nucs[col].pop_cnt  += s[j]->nucs[col].num_indv > 0 ? 1 : 0;
        }

        for (uint j = 0; j < this->_pop_cnt; j++) {
            //
            // Sum the most frequent allele across populations.
            //
            if (s[j]->nucs[col].p_nuc == mp->nucs[col].p_allele)
                mp->nucs[col].p_freq +=
                    s[j]->nucs[col].p * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
            else
                mp->nucs[col].p_freq +=
                    (1 - s[j]->nucs[col].p) * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
            //
            // Sum observed heterozygosity across populations.
            //
            mp->nucs[col].obs_het +=
                s[j]->nucs[col].obs_het * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
        }

        //
        // We want to report the most frequent allele as the P allele. Reorder the alleles
        // if necessary.
        // XXX Possibly unstable for p_freq ~ 0.5. @Nick (July 2016)
        //
        if (mp->nucs[col].p_freq < 0.5) {
            char a = mp->nucs[col].p_allele;
            mp->nucs[col].p_allele = mp->nucs[col].q_allele;
            mp->nucs[col].q_allele = a;
            mp->nucs[col].p_freq   = 1 - mp->nucs[col].p_freq;
            uint b = p_cnt;
            p_cnt = q_cnt;
            q_cnt = b;
        }

        //
        // Check if this is a private allele. Either the site is variable and
        // the allele exists in one population, or the site is fixed and one
        // population is homozygous for the private allele.
        //
        variable_pop = -1;

        if (p_cnt == 1 && q_cnt > 1) {
            for (uint j = 0; j < this->_pop_cnt; j++)
                if (s[j]->nucs[col].p_nuc == mp->nucs[col].p_allele ||
                    s[j]->nucs[col].q_nuc == mp->nucs[col].p_allele)
                    variable_pop = j;
        } else if (p_cnt > 1 && q_cnt == 1) {
            for (uint j = 0; j < this->_pop_cnt; j++)
                if (s[j]->nucs[col].p_nuc == mp->nucs[col].q_allele ||
                    s[j]->nucs[col].q_nuc == mp->nucs[col].q_allele)
                    variable_pop = j;
        }
        mp->nucs[col].priv_allele = variable_pop;
    }

    return 0;
}

int
LocPopSum::tally_ref_alleles(int snp_index, uint16_t &allele_cnt,
                             char &p_allele, char &q_allele,
                             uint16_t &p_cnt, uint16_t &q_cnt)
{
    int  nucs[4] = {0};
    char nuc[2];

    p_allele   = 0;
    q_allele   = 0;
    allele_cnt = 0;

    for (uint j = 0; j < this->_pop_cnt; j++) {
        nuc[0] = 0;
        nuc[1] = 0;
        nuc[0] = this->_per_pop[j]->nucs[snp_index].p_nuc;
        nuc[1] = this->_per_pop[j]->nucs[snp_index].q_nuc;

        for (uint k = 0; k < 2; k++)
            switch(nuc[k]) {
            case 'A':
            case 'a':
                nucs[0]++;
                break;
            case 'C':
            case 'c':
                nucs[1]++;
                break;
            case 'G':
            case 'g':
                nucs[2]++;
                break;
            case 'T':
            case 't':
                nucs[3]++;
                break;
            }
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
        p_allele = 0;
        q_allele = 0;
        return 0;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    i = 0;
    while (p_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 0:
                p_allele = 'A';
                break;
            case 1:
                p_allele = 'C';
                break;
            case 2:
                p_allele = 'G';
                break;
            case 3:
                p_allele = 'T';
                break;
            }
        }
        i++;
    }
    while (q_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 1:
                q_allele = 'C';
                break;
            case 2:
                q_allele = 'G';
                break;
            case 3:
                q_allele = 'T';
                break;
            }
        }
        i++;
    }

    //
    // Tabulate the number of populations the p_allele and the q_allele occur in.
    //
    p_cnt = 0;
    q_cnt = 0;

    for (uint j = 0; j < this->_pop_cnt; j++) {
        nuc[0] = 0;
        nuc[1] = 0;
        nuc[0] = this->_per_pop[j]->nucs[snp_index].p_nuc;
        nuc[1] = this->_per_pop[j]->nucs[snp_index].q_nuc;

        for (uint k = 0; k < 2; k++)
            if (nuc[k] != 0 && nuc[k] == p_allele)
                p_cnt++;
            else if (nuc[k] != 0 && nuc[k] == q_allele)
                q_cnt++;
    }

    return 1;
}

int
LocPopSum::tally_observed_haplotypes(const vector<char *> &obshap, int snp_index)
{
    int  nucs[4] = {0};
    char nuc;

    //
    // Pull each allele for this SNP from the observed haplotype.
    //
    for (uint j = 0; j < obshap.size(); j++) {
        nuc = obshap[j][snp_index];

        switch(nuc) {
        case 'A':
        case 'a':
            nucs[0]++;
            break;
        case 'C':
        case 'c':
            nucs[1]++;
            break;
        case 'G':
        case 'g':
            nucs[2]++;
            break;
        case 'T':
        case 't':
            nucs[3]++;
            break;
        }
    }

    int allele_cnt = 0;
    for (int i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    return allele_cnt;
}

int
LocPopSum::calc_hapstats(const CSLocus *cloc, const Datum **d, const MetaPopInfo &mpopi)
{
    const vector<Pop>  &pops = mpopi.pops();

    for (uint j = 0; j < pops.size(); j++) {

        if ( (this->_hapstats_per_pop[j] = this->haplotype_diversity(pops[j].first_sample, pops[j].last_sample, d)) != NULL) {
            this->_hapstats_per_pop[j]->loc_id = cloc->id;
            this->_hapstats_per_pop[j]->bp     = cloc->sort_bp();

            if (calc_hwp && this->_hapstats_per_pop[j]->hap_cnt > 1) {
                //
                // Initialize the object and Calculate Hardy-Weinberg Proportions.
                //
                GuoThompson_Hwp *hwp = new GuoThompson_Hwp(this->_hapstats_per_pop[j]->hap_cnt);

                hwp->exec_locus(pops[j].first_sample, pops[j].last_sample, d, this->_hapstats_per_pop[j]->hap_cnt);

                this->_hapstats_per_pop[j]->stat[2] = hwp->_p_value;
                this->_hapstats_per_pop[j]->stat[3] = hwp->_se;

                delete hwp;
            }
        }
    }

    return 0;
}

LocStat *
LocPopSum::haplotype_diversity(int start, int end, Datum const*const* d)
{
    map<string, double>::iterator hit;
    vector<string>      haplotypes;
    map<string, double> hap_freq;
    map<string, int>    hap_index;
    double gene_diversity = 0.0;
    double hapl_diversity = 0.0;

    LocStat *lstat = new LocStat;

    //
    // Tabulate the haplotypes in this population.
    //
    double n = count_haplotypes_at_locus(start, end, d, hap_freq);

    //
    // If this haplotype is fixed, don't calculate any statistics.
    //
    if (n == 0)
        return NULL;

    //
    // Store a summary of the haplotype counts to output below.
    //
    stringstream sstr;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++)
        sstr << hit->first << ":" << hit->second  << ";";
    lstat->hap_str = sstr.str().substr(0, sstr.str().length() - 1);

    //
    // Determine an ordering for the haplotypes. Convert haplotype counts into frequencies.
    //
    uint k = 0;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++) {
        hap_index[hit->first] = k;
        haplotypes.push_back(hit->first);
        k++;

        // cout << "  Haplotype '" << hit->first << "' occured " << hit->second  << " times; ";

        hit->second = hit->second / n;

        // cout << " frequency of " << hit->second << "%\n";
    }
    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists = new double *[hap_index.size()];
    for (k = 0; k < hap_index.size(); k++) {
        hdists[k] = new double[hap_index.size()];
        memset(hdists[k], 0, hap_index.size() * sizeof(double));
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(hap_index, hdists);

    //
    // Calculate haplotype diversity, Pi.
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
        for (uint j = 0; j < haplotypes.size(); j++) {
            hapl_diversity +=
                hap_freq[haplotypes[i]] *
                hap_freq[haplotypes[j]] *
                hdists[hap_index[haplotypes[i]]][hap_index[haplotypes[j]]];
        }
    }
    hapl_diversity = (n / (n-1)) * hapl_diversity;

    //
    // Calculate gene diversity.
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
        gene_diversity += hap_freq[haplotypes[i]] * hap_freq[haplotypes[i]];
    }
    gene_diversity = (n / (n - 1)) * (1 - gene_diversity);

    lstat->alleles = n;
    lstat->stat[0] = gene_diversity;
    lstat->stat[1] = hapl_diversity;
    lstat->hap_cnt = haplotypes.size();

    // cout << "  Population " << pop_id << " has haplotype diversity (pi) of " << s[pop_index]->pi << "\n";

    for (k = 0; k < hap_index.size(); k++)
        delete [] hdists[k];
    delete [] hdists;

    return lstat;
}

int
nuc_substitution_dist(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    const char *p, *q;
    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {

            dist = 0.0;
            p    = haplotypes[i].c_str();
            q    = haplotypes[j].c_str();

            while (*p != '\0' && *q != '\0') {
                if (*p != *q) dist++;
                p++;
                q++;
            }

            hdists[i][j] = dist;
            hdists[j][i] = dist;
        }
    }

    // //
    // // Print the distance matrix.
    // //
    // cout << "  ";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++)
    //  cout << "\t" << hit->first;
    // cout << "\n";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
    //  cout << "  " << hit->first;
    //  for (hit_2 = loc_hap_index.begin(); hit_2 != loc_hap_index.end(); hit_2++)
    //      cout << "\t" << hdists[hit->second][hit_2->second];
    //  cout << "\n";
    // }
    // cout << "\n";

    return 0;
}

bool
uncalled_haplotype(const char *haplotype)
{
    for (const char *p = haplotype; *p != '\0'; p++)
        if (*p == 'N' || *p == 'n')
            return true;
    return false;
}

double
count_haplotypes_at_locus(int start, int end, Datum const*const* d, map<string, double> &hap_cnts)
{
    double n = 0.0;

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
            if(!uncalled_haplotype(d[i]->obshap[0])) {
                n += 2;
                hap_cnts[d[i]->obshap[0]] += 2;
            }
        } else {
            // Heterozygote.
            for (uint j = 0; j < d[i]->obshap.size(); j++) {
                if(!uncalled_haplotype(d[i]->obshap[j])) {
                    n++;
                    hap_cnts[d[i]->obshap[j]]++;
                }
            }
        }
    }
    return n;
}

LocusDivergence::LocusDivergence(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++)
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            this->_mean_fst.push_back(0.0);
            this->_mean_fst_cnt.push_back(0.0);
            this->_mean_phist.push_back(0.0);
            this->_mean_phist_cnt.push_back(0.0);
            this->_mean_fstp.push_back(0.0);
            this->_mean_fstp_cnt.push_back(0.0);
            this->_mean_dxy.push_back(0.0);
            this->_mean_dxy_cnt.push_back(0.0);
        }
}

int
LocusDivergence::clear(const vector<LocBin *> &loci)
{
    for (uint i = 0; i < this->_snps.size(); i++) {
        assert (this->_snps[i].size() == loci.size());
        for (uint j = 0; j < this->_snps[i].size(); j++) {
            for (uint k = 0; k < loci[j]->cloc->len; k++)
                if (this->_snps[i][j][k] != NULL)
                    delete this->_snps[i][j][k];
            delete [] this->_snps[i][j];
        }
        this->_snps[i].clear();
    }
    this->_snps.clear();

    for (uint i = 0; i < this->_haplotypes.size(); i++) {
        for (uint j = 0; j < this->_haplotypes[i].size(); j++)
            delete this->_haplotypes[i][j];
        this->_haplotypes[i].clear();
    }
    this->_haplotypes.clear();

    for (uint i = 0; i < this->_metapop_haplotypes.size(); i++)
        if (this->_metapop_haplotypes[i] != NULL)
            delete this->_metapop_haplotypes[i];
    this->_metapop_haplotypes.clear();

    return 0;
}

int
LocusDivergence::snp_divergence(const vector<LocBin *> &loci)
{
    uint i = 0;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {

        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {

            double fmean = 0.0;
            double fcnt  = 0.0;

            //
            // Preallocate our vector of PopPairs so each thread can return its results in the proper order.
            //
            this->_snps.push_back(vector<PopPair **>(loci.size(), NULL));

            #pragma omp parallel
            {
                vector<PopPair **> *pairs = &this->_snps.back();

                #pragma omp for schedule(dynamic, 1) reduction(+:fmean, fcnt)
                for (uint k = 0; k < loci.size(); k++) {

                    const LocBin *loc = (const LocBin *) loci[k];
                    uint     cloc_len = loc->cloc->len;
                    PopPair      **pp = new PopPair *[cloc_len];

                    for (uint pos = 0; pos < cloc_len; pos++) {
                        pp[pos] = this->Fst(loc->cloc, loc->s, pop_1, pop_2, pos);

                        //
                        // Locus is fixed in both populations, or was only found in one population.
                        //
                        if (pp[pos]->pi == 0) {
                            delete pp[pos];
                            pp[pos] = NULL;
                            continue;
                        }

                        pp[pos]->loc_id = loc->cloc->id;
                        pp[pos]->bp     = loc->cloc->sort_bp(pos);
                        pp[pos]->col    = pos;
                        pp[pos]->pop_1  = pop_1;
                        pp[pos]->pop_2  = pop_2;

                        //
                        // Apply user-selected correction to the Fst values.
                        //
                        switch(fst_correction) {
                        case p_value:
                            if (pp[pos] != NULL) {
                                pp[pos]->stat[0] = pp[pos]->fet_p < p_value_cutoff ? pp[pos]->fst : 0;
                                pp[pos]->stat[1] = pp[pos]->fet_p < p_value_cutoff ? pp[pos]->amova_fst : 0;
                            }
                            break;
                        case no_correction:
                        case bonferroni_win:
                        case bonferroni_gen:
                            if (pp[pos] != NULL) {
                                pp[pos]->stat[0] = pp[pos]->fst;
                                pp[pos]->stat[1] = pp[pos]->amova_fst;
                            }
                            break;
                        }

                        fmean += pp[pos]->amova_fst;
                        fcnt++;
                    }

                    pairs->at(k) = pp;
                }
            }

            this->_mean_fst[i]     += fmean;
            this->_mean_fst_cnt[i] += fcnt;
            i++;
        }
    }

    return 0;
}

PopPair *
LocusDivergence::Fst(const CSLocus *cloc, const LocPopSum *s, int pop_1, int pop_2, int pos)
{
    const LocSum *s_1 = s->per_pop(pop_1);
    const LocSum *s_2 = s->per_pop(pop_2);

    PopPair *pair = new PopPair();

    //
    // If this locus only appears in one population do not calculate Fst.
    //
    if (s_1->nucs[pos].num_indv == 0 || s_2->nucs[pos].num_indv == 0)
        return pair;

    //
    // Calculate Fst at a locus, sub-population relative to that found in the entire population
    //   Fst = 1 - (Sum_j( (n_j choose 2) * pi_j)) / (pi_all * Sum_j( (n_j choose 2) ))
    //
    double n_1, n_2, pi_1, pi_2;

    n_1  = s_1->nucs[pos].num_indv * 2;
    n_2  = s_2->nucs[pos].num_indv * 2;
    pi_1 = s_1->nucs[pos].pi;
    pi_2 = s_2->nucs[pos].pi;

    if (pi_1 == 0 && pi_2 == 0 && s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc)
        return pair;

    //
    // Calculate Pi over the entire pooled population.
    //
    // First, make sure this site is compatible between the two populations (no more than two alleles).
    //
    char nucs[4];
    int  ncnt[4] = {0};

    nucs[0] = s_1->nucs[pos].p_nuc;
    nucs[1] = s_1->nucs[pos].q_nuc;
    nucs[2] = s_2->nucs[pos].p_nuc;
    nucs[3] = s_2->nucs[pos].q_nuc;

    for (int i = 0; i < 4; i++)
        switch(nucs[i]) {
        case 'A':
            ncnt[0]++;
            break;
        case 'C':
            ncnt[1]++;
            break;
        case 'G':
            ncnt[2]++;
            break;
        case 'T':
            ncnt[3]++;
            break;
        }

    int allele_cnt = 0;
    for (int i = 0; i < 4; i++)
        if (ncnt[i] > 0) allele_cnt++;

    if (allele_cnt > 2)
        return pair;

    double tot_alleles = n_1 + n_2;
    double p_1 = round(n_1 * s_1->nucs[pos].p);
    double q_1 = n_1 - p_1;
    double p_2 =
        s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc ?
        s_2->nucs[pos].p : (1 - s_2->nucs[pos].p);
    p_2 = round(n_2 * p_2);
    double q_2 = n_2 - p_2;

    double pi_all = s->pi(tot_alleles, p_1 + p_2, q_1 + q_2);

    double bcoeff_1 = s->binomial_coeff(n_1, 2);
    double bcoeff_2 = s->binomial_coeff(n_2, 2);

    double num = (bcoeff_1 * pi_1) + (bcoeff_2 * pi_2);
    double den = pi_all * (bcoeff_1 + bcoeff_2);

    double Fst = 1 - (num / den);

    pair->alleles = tot_alleles;
    pair->fst     = Fst;
    pair->pi      = pi_all;

    this->fishers_exact_test(pair, p_1, q_1, p_2, q_2);

    // cout << "Locus: " << locus << ", pos: " << pos << "\n"
    //   << "    p_1.nuc: " << s_1->nucs[pos].p_nuc << "; q_1.nuc: " << s_1->nucs[pos].q_nuc
    //   << "; p_2.nuc: " << s_2->nucs[pos].p_nuc << "; q_2.nuc: " << s_2->nucs[pos].q_nuc << "\n"
    //   << "    Total alleles: " << tot_alleles << "; " << " s_1.p: " << s_1->nucs[pos].p
    //   << "; s_2.p: " << s_2->nucs[pos].p << "\n"
    //   << "    p_1: " << p_1 << "; q_1: " << q_1 << " p_2: " << p_2 << "; q_2: " << q_2 << "\n"
    //   << "    Pi1: " << pi_1 << "; Pi2: " << pi_2 << "; PiAll: " << pi_all << "\n"
    //   << "    N1: " << n_1 << "; N1 choose 2: " << bcoeff_1 << "\n"
    //   << "    N2: " << n_2 << "; N2 choose 2: " << bcoeff_2 << "\n"
    //   << "  Fst: " << Fst << "\n";

    //
    // Calculate Fst (corrected for different samples sizes) using an AMOVA method,
    // correcting for unequal sample sizes.
    // Derived from Weir, _Genetic Data Analysis II_, chapter 5, "F Statistics,", pp166-167.
    //
    double p_1_freq = s_1->nucs[pos].p;
    double q_1_freq = 1 - p_1_freq;
    double p_2_freq =
        s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc ?
        s_2->nucs[pos].p : (1 - s_2->nucs[pos].p);
    double q_2_freq = 1 - p_2_freq;

    double p_avg_cor =
        ( (s_1->nucs[pos].num_indv * p_1_freq) + (s_2->nucs[pos].num_indv * p_2_freq) ) /
        ( s_1->nucs[pos].num_indv + s_2->nucs[pos].num_indv );
    double n_avg_cor = 2 * ((s_1->nucs[pos].num_indv / 2) + (s_2->nucs[pos].num_indv / 2));

    pair->amova_fst =
        (
         (s_1->nucs[pos].num_indv * pow((p_1_freq - p_avg_cor), 2) +
          s_2->nucs[pos].num_indv * pow((p_2_freq - p_avg_cor), 2))
         /
         n_avg_cor
         )
        /
        (p_avg_cor * (1 - p_avg_cor));

    if (log_fst_comp) {
        pair->comp     = new double[18];
        pair->comp[0]  = n_1;
        pair->comp[1]  = n_2;
        pair->comp[2]  = tot_alleles;
        pair->comp[3]  = p_1;
        pair->comp[4]  = q_1;
        pair->comp[5]  = p_2;
        pair->comp[6]  = q_2;
        pair->comp[7]  = pi_1;
        pair->comp[8]  = pi_2;
        pair->comp[9]  = pi_all;
        pair->comp[10] = bcoeff_1;
        pair->comp[11] = bcoeff_2;
        pair->comp[12] = p_1_freq;
        pair->comp[13] = q_1_freq;
        pair->comp[14] = p_2_freq;
        pair->comp[15] = q_2_freq;
        pair->comp[16] = p_avg_cor;
        pair->comp[17] = n_avg_cor;
    }

    // //
    // // Calculate Fst using a pure parametric method (assumes allele counts are real, not
    // // samples). Jakobsson, Edge, and Rosenberg. "The Relationship Between Fst and the
    // // Frequency of the Most Frequent Allele." Genetics 193:515-528. Equation 4.
    // //
    // double sigma_1 = p_1_freq + q_1_freq;
    // double sigma_2 = p_2_freq + q_2_freq;
    // double delta_1 = fabs(p_1_freq - p_2_freq);
    // double delta_2 = fabs(q_1_freq - q_2_freq);

    // pair->jakob_fst = (pow(delta_1, 2) + pow(delta_2, 2)) / ( 4 - (pow(sigma_1, 2) + pow(sigma_2, 2)) );

    return pair;
}

int
LocusDivergence::fishers_exact_test(PopPair *pair, double p_1, double q_1, double p_2, double q_2)
{
    //                            | Allele1 | Allele2 |
    // Fisher's Exact Test:  -----+---------+---------+
    //                       Pop1 |   p_1   |   q_1   |
    //                       Pop2 |   p_2   |   q_2   |
    // Probability p:
    // p = ((p_1 + q_1)!(p_2 + q_2)!(p_1 + p_2)!(q_1 + q_2)!) / (n! p_1! q_1! p_2! q_2!)
    //
    // According to:
    //   Jerrold H. Zar, "A fast and efficient algorithm for the Fisher exact test."
    //   Behavior Research Methods, Instruments, & Computers 1987, 19(4): 43-44
    //
    // Probability p can be calculated as three binomial coefficients:
    //   Let p_1 + q_1 = r_1; p_2 + q_2 = r_2; p_1 + p_2 = c_1; q_1 + q_2 = c_2
    //
    //   p = (r_1 choose p_1)(r_2 choose p_2) / (n choose c_1)
    //
    // Fisher's Exact test algorithm implemented according to Sokal and Rohlf, _Biometry_, section 17.4.
    //

    double r_1  = p_1 + q_1;
    double r_2  = p_2 + q_2;
    double c_1  = p_1 + p_2;
    double d_1  = p_1 * q_2;
    double d_2  = p_2 * q_1;
    double n    = r_1 + r_2;
    double p    = 0.0;

    // char p1_str[32], p2_str[32], q1_str[32], q2_str[32];
    // sprintf(p1_str, "% 3.0f", p_1);
    // sprintf(q1_str, "% 3.0f", q_1);
    // sprintf(p2_str, "% 3.0f", p_2);
    // sprintf(q2_str, "% 3.0f", q_2);
    //
    // cout
    //  << "     | Allele1 | Allele2 | " << "\n"
    //  << "-----+---------+---------+" << "\n"
    //  << "Pop1 |   " << p1_str << "   |   " << q1_str << "   |" << "\n"
    //  << "Pop2 |   " << p2_str << "   |   " << q2_str << "   |" << "\n\n";

    //
    // Compute the first tail.
    //
    double p1     = p_1;
    double q1     = q_1;
    double p2     = p_2;
    double q2     = q_2;
    double tail_1 = 0.0;
    double den    = LocPopSum::binomial_coeff(n, c_1);

    //
    // If (p_1*q_2 - p_2*q_1) < 0 decrease cells p_1 and q_2 by one and add one to p_2 and q_1.
    // Compute p and repeat until one or more cells equal 0.
    //
    if (d_1 - d_2 < 0) {
        do {
            p = (LocPopSum::binomial_coeff(r_1, p1) * LocPopSum::binomial_coeff(r_2, p2)) / den;

            tail_1 += p;
            p1--;
            q2--;
            p2++;
            q1++;
        } while (p1 >= 0 && q2 >= 0);

    } else {
        //
        // Else, if (p_1*q_2 - p_2*q_1) > 0 decrease cells p_2 and q_1 by one and add one to p_1 and q_2.
        // Compute p and repeat until one or more cells equal 0.
        //
        do {
            p = (LocPopSum::binomial_coeff(r_1, p1) * LocPopSum::binomial_coeff(r_2, p2)) / den;

            tail_1 += p;

            p2--;
            q1--;
            p1++;
            q2++;
        } while (p2 >= 0 && q1 >= 0);
    }

    //
    // Compute the second tail.
    //
    double tail_2 = 0.0;
    p = 0;

    //
    // If (p_1*q_2 - p_2*q_1) < 0, set to zero the smaller of the two frequencies, adjusting the other values
    // to keep the marginals the same.
    //
    if (d_1 - d_2 < 0) {
        if (p2 < q1) {
            q2 += p2;
            p1 += p2;
            q1 -= p2;
            p2  = 0;
        } else {
            p1 += q1;
            q2 += q1;
            p2 -= q1;
            q1  = 0;
        }
    } else {
        if (p1 < q2) {
            q1 += p1;
            p2 += p1;
            q2 -= p1;
            p1  = 0;
        } else {
            p2 += q2;
            q1 += q2;
            p1 -= q2;
            q2  = 0;
        }
    }

    //
    // If (p_1*q_2 - p_2*q_1) < 0 decrease cells p_1 and q_2 by one and add one to p_2 and q_1.
    // Compute p and repeat until tail_2 > tail_1.
    //
    if (d_1 - d_2 < 0) {
        do {
            p = (LocPopSum::binomial_coeff(r_1, p1) * LocPopSum::binomial_coeff(r_2, p2)) / den;

            tail_2 += p;

            p1--;
            q2--;
            p2++;
            q1++;
        } while (tail_2 < tail_1 && p1 >= 0 && q2 >= 0);

        tail_2 -= p;

    } else {
        //
        // Else, if (p_1*q_2 - p_2*q_1) > 0 decrease cells p_2 and q_1 by one and add one to p_1 and q_2.
        // Compute p and repeat until one or more cells equal 0.
        //
        do {
            p = (LocPopSum::binomial_coeff(r_1, p1) * LocPopSum::binomial_coeff(r_2, p2)) / den;

            tail_2 += p;

            p2--;
            q1--;
            p1++;
            q2++;
        } while (tail_2 < tail_1 && p2 >= 0 && q1 >= 0);

        tail_2 -= p;
    }

    pair->fet_p = tail_1 + tail_2;

    if (pair->fet_p > 1.0) pair->fet_p = 1.0;

    //
    // Calculate the odds ratio. To account for possible cases were one allele frequency is
    // zero, we will increment all allele frequencies by one.
    //
    if (p_1 == 0 || q_1 == 0 || p_2 == 0 || q_2 == 0) {
        p_1++;
        q_1++;
        p_2++;
        q_2++;
    }
    pair->fet_or = (p_1 * q_2) / (q_1 * p_2);

    double ln_fet_or = pair->fet_or > 0 ? log(pair->fet_or) : 0.0;

    //
    // Calculate the standard error of the natural log of the odds ratio
    //
    double se = pair->fet_or > 0 ? sqrt((1 / p_1) + (1 / q_1) + (1 / p_2) + (1 / q_2)) : 0.0;

    //
    // Calculate the confidence intervals of the natural log of the odds ratio.
    //
    double ln_ci_low  = pair->fet_or > 0 ? ln_fet_or - (1.96 * se) : 0;
    double ln_ci_high = pair->fet_or > 0 ? ln_fet_or + (1.96 * se) : 0;

    //
    // Convert the confidence intervals out of natural log space
    //
    pair->ci_low  = pair->fet_or > 0 ? exp(ln_ci_low)  : 0;
    pair->ci_high = pair->fet_or > 0 ? exp(ln_ci_high) : 0;
    pair->lod     = fabs(log10(pair->fet_or));

    return 0;
}

int
LocusDivergence::haplotype_divergence_pairwise(const vector<LocBin *> &loci)
{
    uint i = 0;
    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {

            vector<int> subpop_ids;
            subpop_ids.push_back(pop_1);
            subpop_ids.push_back(pop_2);

            this->_haplotypes.push_back(vector<HapStat *>(loci.size(), NULL));

            double pmean   = 0.0;
            double fmean   = 0.0;
            double dxymean = 0.0;
            double cnt     = 0.0;
                
            #pragma omp parallel
            {
                vector<HapStat *> *haps = &this->_haplotypes.back();

                #pragma omp for schedule(dynamic, 1) reduction(+:pmean, fmean, cnt)
                for (uint k = 0; k < loci.size(); k++) {
                    const CSLocus *loc = (const CSLocus *) loci[k]->cloc;
                    const LocSum **s   = loci[k]->s->all_pops();
                    const Datum  **d   = (const Datum **) loci[k]->d;

                    HapStat *h;
                    //
                    // If this locus only appears in one population or there is only a single haplotype,
                    // do not calculate haplotype F stats.
                    //
                    if (loc->snps.size() == 0 || fixed_locus(d, subpop_ids))
                        h = NULL;
                    else
                        h = this->haplotype_amova(d, s, subpop_ids);

                    if (h != NULL) {
                        h->pop_1   = pop_1;
                        h->pop_2   = pop_2;
                        h->stat[4] = haplotype_d_est(d, s, subpop_ids);
                        h->stat[5] = haplotype_dxy(d, loc->len, subpop_ids);

                        h->loc_id = loc->id;
                        h->bp     = loc->sort_bp();

                        pmean   += h->stat[0];
                        fmean   += h->stat[3];
                        dxymean += h->stat[5];
                        cnt++;
                    }

                    haps->at(k) = h;
                }
            }

            this->_mean_phist[i]     += pmean;
            this->_mean_phist_cnt[i] += cnt;
            this->_mean_fstp[i]      += fmean;
            this->_mean_fstp_cnt[i]  += cnt;
            this->_mean_dxy[i]       += dxymean;
            this->_mean_dxy_cnt[i]   += cnt;
            i++;
        }
    }

    return 0;
}

int
LocusDivergence::haplotype_divergence(const vector<LocBin *> &loci)
{
    //
    // Create a list of all the populations we have.
    //
    vector<int> pop_ids;
    for (size_t i = 0; i < this->_mpopi->pops().size(); ++i)
        pop_ids.push_back(i);

    this->_metapop_haplotypes.resize(loci.size(), NULL);

    #pragma omp parallel
    {
        vector<HapStat *> *haps = &this->_metapop_haplotypes;

        #pragma omp for schedule(dynamic, 1)
        for (uint k = 0; k < loci.size(); k++) {
            const CSLocus *loc = (const CSLocus *) loci[k]->cloc;
            const LocSum **s   = loci[k]->s->all_pops();
            const Datum  **d   = (const Datum **) loci[k]->d;

            HapStat *h;
            //
            // If this locus only appears in one population or there is only a single haplotype,
            // do not calculate haplotype F stats.
            //
            if (loc->snps.size() == 0 || fixed_locus(d, pop_ids))
                h = NULL;
            else
                h = this->haplotype_amova(d, s, pop_ids);

            if (h != NULL) {
                h->stat[4] = this->haplotype_d_est(d, s, pop_ids);

                h->loc_id  = loc->id;
                h->bp      = loc->sort_bp();
            }

            haps->at(k) = h;
        }
    }

    return 0;
}

HapStat *
LocusDivergence::haplotype_amova(const Datum **d, const LocSum **s, vector<int> &pop_ids)
{
    map<string, int>          loc_hap_index;
    vector<string>            loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;
    map<int, vector<int> >    grp_members;
    vector<int>               grps;

    map<string, int>::iterator hit, hit_2;

    HapStat  *h;

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int pop_id : pop_ids) {
        const Pop& pop = this->_mpopi->pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) continue;

            if (d[i]->obshap.size() > 2) {
                continue;

            } else if (d[i]->obshap.size() == 1) {
                if(!uncalled_haplotype(d[i]->obshap[0])) {
                    loc_hap_index[d[i]->obshap[0]]++;
                    loc_haplotypes.push_back(d[i]->obshap[0]);
                    loc_haplotypes.push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                }
            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    if(!uncalled_haplotype(d[i]->obshap[0])) {
                        loc_hap_index[d[i]->obshap[j]]++;
                        loc_haplotypes.push_back(d[i]->obshap[j]);
                        pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
                    }
                }
            }
        }
    }

    //
    // What is the total number of populations that had valid haplotypes.
    //
    double valid_pop_cnt = 0.0;
    for (int pop_id : pop_ids) {
        if (pop_haplotypes[pop_id].size() > 0)
            valid_pop_cnt++;
    }

    //
    // If we filtered a population out at this locus make sure that we still have at least one
    // representative present in each group.
    //
    set<int> uniq_grps;
    for (size_t pop_id = 0; pop_id < this->_mpopi->pops().size(); ++pop_id) {
        const Pop& pop = this->_mpopi->pops()[pop_id];
        if (pop_haplotypes.count(pop_id) > 0) {
            uniq_grps.insert(pop.group);
            grp_members[pop.group].push_back(pop_id);
        }
    }
    set<int>::iterator uit;
    for (uit = uniq_grps.begin(); uit != uniq_grps.end(); uit++)
        grps.push_back(*uit);

    if (grps.size() == 0)
        return NULL;

    // cout << "Groups: ";
    // for (uint i = 0; i < grps.size(); i++)
    //     cout << grps[i] << ", ";
    // cout << "\n";
    // for (git = grp_members.begin(); git != grp_members.end(); git++) {
    //     cout << "Group " << git->first << ": ";
    //     for (uint i = 0; i < git->second.size(); i++)
    //  cout << git->second[i] << ", ";
    //     cout << "\n";
    // }

    //
    // Determine an ordering for the haplotypes.
    //
    uint m = 0;
    for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
        loc_hap_index[hit->first] = m;
        m++;
    }

    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists     = new double *[loc_hap_index.size()];
    double **hdists_max = new double *[loc_hap_index.size()];
    for (uint k = 0; k < loc_hap_index.size(); k++) {
        hdists[k] = new double[loc_hap_index.size()];
        memset(hdists[k], 0, loc_hap_index.size() * sizeof(double));
        hdists_max[k] = new double[loc_hap_index.size()];
        memset(hdists_max[k], 0, loc_hap_index.size() * sizeof(double));
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(loc_hap_index, hdists);

    //
    // Calculate the sum of squared distances in each subset: total, within populations, across populations
    // and withing groups, and across groups.
    //
    double ssd_total = this->amova_ssd_total(loc_haplotypes, loc_hap_index, hdists);
    double ssd_wp    = this->amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    double ssd_ap_wg = this->amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    double ssd_ag    = grps.size() > 1 ? this->amova_ssd_ag(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, ssd_total) : 0.0;

    //
    // Calculate n
    //
    double n        = 0.0;
    double n_1      = 0.0;
    double n_2      = 0.0;
    double s_g      = 0.0;
    double tot_cnt  = 0.0;
    double grp_cnt  = 0.0;
    double num_grps = grps.size();
    double a        = 0.0;
    double b        = 0.0;

    for (uint g = 0; g < num_grps; g++) {
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            tot_cnt += (double) pop_haplotypes[pop_id_1].size();
        }
    }
    for (uint g = 0; g < num_grps; g++) {
        grp_cnt = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            grp_cnt += (double) pop_haplotypes[pop_id_1].size();
        }

        a = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            a += (double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / grp_cnt;
        }
        s_g += a;
    }
    n = (tot_cnt - s_g) / (double) (valid_pop_cnt - num_grps);

    // cout << "  n: "<< n << "\n";

    if (num_grps > 1) {
        //
        // Calculate n'
        //
        a = 0.0;
        for (uint g = 0; g < num_grps; g++) {
            for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
                int pop_id_1 = grp_members[grps[g]][r];
                a += ((double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / tot_cnt);
            }
        }
        n_1 = (s_g - a) / (double) (num_grps - 1.0);

        // cout << "  n': "<< n_1 << "\n";

        //
        // Calculate n''
        //
        for (uint g = 0; g < num_grps; g++) {
            a = 0.0;
            for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
                int pop_id_1 = grp_members[grps[g]][r];
                a += pop_haplotypes[pop_id_1].size();
            }
            b += ((a * a) / tot_cnt);
        }
        n_2 = (tot_cnt - b) / (double) (num_grps - 1);

        // cout << "  n'': "<< n_2 << "\n";
    }

    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    double msd_ag    = num_grps > 1 ? ssd_ag / (double) (num_grps - 1) : 0.0;
    double msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    double msd_wp    = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    double msd_total = ssd_total / ((double) (loc_haplotypes.size() - 1));

    double sigma_c     = msd_wp;
    double sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    double sigma_a     = 0.0;

    if (grps.size() > 1)
        sigma_a = (msd_ag - sigma_c - (n_1 * sigma_b)) / n_2;

    // Arlequin seems to sum the variance components instead of independently calculating sigma_total: MSD(total) = SSD(total)/degrees.of.freedom
    double sigma_total = sigma_a + sigma_b + sigma_c; // msd_total;

    double phi_st = 0.0;
    double phi_ct = 0.0;
    double phi_sc = 0.0;

    if (grps.size() > 1) {
        phi_st = sigma_total > 0.0 ? (sigma_a + sigma_b) / sigma_total : 0.0;
        phi_ct = sigma_total > 0.0 ?  sigma_a / sigma_total : 0.0;
        phi_sc = (sigma_a + sigma_b) > 0.0 ?  sigma_b / (sigma_b + sigma_c) : 0.0;
    } else {
        phi_st = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    }

    // cout << "  MSD(AG): " << msd_ag  << "; MSD(AP/WG): " << msd_ap_wg << "; MSD(WP): " << msd_wp  << "; MSD(TOTAL): "  << msd_total   << "\n"
    //      << "  Sigma_a: " << sigma_a << "; Sigma_b: "    << sigma_b   << "; Sigma_c: " << sigma_c << "; Sigma_Total: " << sigma_total << "\n"
    //      << "  Phi_st: "  << phi_st  << "; Phi_ct: "     << phi_ct    << "; Phi_sc: "  << phi_sc  << "\n";

    //
    // Calculate Fst' = Fst / Fst_max
    //
    // First calculate Fst.
    //
    // To calculate Fst instead of Phi_st, we need to reset our distance matrix to return 1 if haplotypes are different, 0 otherwise.
    //
    nuc_substitution_identity(loc_hap_index, hdists);
    ssd_wp    = amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    msd_ap_wg   = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    msd_wp      = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    sigma_c     = msd_wp;
    sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    sigma_total = sigma_b + sigma_c;

    double fst = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;

    //
    // Now calculate Fst_max.
    //
    // Reset our distance matrix to give maximum possible distance between haplotypes
    // and recalculate sum of squared deviations across groups.
    //
    nuc_substitution_identity_max(loc_hap_index, hdists_max);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists_max);

    //
    // Recalculate the mean square deviations, given maximum divergence between populations.
    //
    msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    sigma_b   = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;

    double fst_max = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    double fst_1   = fst_max > 0.0     ? fst / fst_max         : 0.0;

    //
    // Cache the results so we can print them in order below, once the parallel code has executed.
    //
    h = new HapStat;
    h->alleles = tot_cnt;
    h->popcnt  = valid_pop_cnt;

    if (log_fst_comp) {
        h->comp = new double[15];
        h->comp[0]  = ssd_wp;
        h->comp[1]  = ssd_ap_wg;
        h->comp[2]  = ssd_ag;
        h->comp[3]  = ssd_total;
        h->comp[4]  = msd_wp;
        h->comp[5]  = msd_ap_wg;
        h->comp[6]  = msd_ag;
        h->comp[7]  = msd_total;
        h->comp[8]  = n;
        h->comp[9]  = n_1;
        h->comp[10] = n_2;
        h->comp[11] = sigma_a;
        h->comp[12] = sigma_b;
        h->comp[13] = sigma_c;
        h->comp[14] = sigma_total;
    }

    h->stat[0] = phi_st;
    h->stat[1] = phi_ct;
    h->stat[2] = phi_sc;
    h->stat[3] = fst_1;

    for (uint k = 0; k < loc_hap_index.size(); k++) {
        delete [] hdists[k];
        delete [] hdists_max[k];
    }
    delete [] hdists;
    delete [] hdists_max;

    return h;
}

double
LocusDivergence::amova_ssd_total(vector<string> &loc_haplotypes, map<string, int> &loc_hap_index, double **hdists)
{
    //
    // Calculate sum of squared deviations for the total sample, SSD(Total)
    //
    double ssd_total = 0.0;

    for (uint j = 0; j < loc_haplotypes.size(); j++) {
        for (uint k = 0; k < loc_haplotypes.size(); k++) {
            ssd_total += hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]];
            // cout << j << "\t"
            //   << k << "\t"
            //   << loc_haplotypes[j] << "\t"
            //   << loc_haplotypes[k] << "\t"
            //   << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
        }
    }
    ssd_total = (1.0 / (double) (2*loc_haplotypes.size())) * ssd_total;
    // cout << "  ssd_total: "<< ssd_total << "\n";

    return ssd_total;
}

double
LocusDivergence::amova_ssd_wp(vector<int> &grps, map<int, vector<int>> &grp_members,
                              map<string, int> &loc_hap_index, map<int, vector<string>> &pop_haplotypes,
                              double **hdists)
{
    //
    // Calculate the sum of squared deviations within populations, SSD(WP)
    //
    double ssd_wp = 0.0;
    double ssd    = 0.0;
    int    pop_id;

    for (uint g = 0; g < grps.size(); g++) {
        for (uint i = 0; i < grp_members[grps[g]].size(); i++) {
            pop_id = grp_members[grps[g]][i];
            ssd = 0.0;

            for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
                for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
                    ssd += hdists[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
                    // cout << pop_id << "\t"
                    //   << j << "\t"
                    //   << k << "\t"
                    //   << loc_haplotypes[j] << "\t"
                    //   << loc_haplotypes[k] << "\t"
                    //   << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
                }
            }

            if (pop_haplotypes[pop_id].size() > 0)
                ssd_wp += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
        }
    }
    // cout << "  ssd_wp: "<< ssd_wp << "\n";

    return ssd_wp;
}

double
LocusDivergence::amova_ssd_ap_wg(vector<int> &grps, map<int, vector<int>> &grp_members,
                                 map<string, int> &loc_hap_index, map<int, vector<string>> &pop_haplotypes,
                                 double **hdists_1, double **hdists_2)
{
    //
    // Calculate the sum of squared deviations across populations and within groups, SSD(AP/WG)
    //
    double ssd_ap_wg = 0.0;
    double ssd       = 0.0;
    double ssd_1     = 0.0;
    double ssd_2     = 0.0;
    double den       = 0.0;
    int    pop_id, pop_id_1, pop_id_2;

    for (uint g = 0; g < grps.size(); g++) {

        ssd_1 = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];

            for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

                for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
                    pop_id_2 = grp_members[grps[g]][s];

                    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
                        if (pop_id_1 == pop_id_2)
                            ssd_1 += hdists_1[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                        else
                            ssd_1 += hdists_2[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                    }
                }
            }
        }

        den = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];
            den += 2 * pop_haplotypes[pop_id_1].size();
        }

        ssd_1 = ssd_1 / den;

        ssd_2 = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id = grp_members[grps[g]][r];
            ssd = 0.0;

            for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
                for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
                    ssd += hdists_1[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
                }
            }

            if (pop_haplotypes[pop_id].size() > 0)
                ssd_2 += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
        }

        ssd_ap_wg += ssd_1 - ssd_2;
    }
    // cout << "  ssd_ap_wg: "<< ssd_ap_wg << "\n";

    return ssd_ap_wg;
}

double
LocusDivergence::amova_ssd_ag(vector<int> &grps, map<int, vector<int>> &grp_members,
                              map<string, int> &loc_hap_index, map<int, vector<string>> &pop_haplotypes,
                              double **hdists, double ssd_total)
{
    //
    // Calculate the sum of squared deviations across groups, SSD(AG)
    //
    int    pop_id_1, pop_id_2;
    double ssd_ag = 0.0;
    double ssd    = 0.0;
    double ssd_1  = 0.0;
    double den    = 0.0;

    for (uint g = 0; g < grps.size(); g++) {
        ssd_1 = 0.0;

        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];

            for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

                for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
                    pop_id_2 = grp_members[grps[g]][s];

                    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
                        ssd_1 += hdists[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                    }
                }
            }
        }

        den = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];
            den += 2 * pop_haplotypes[pop_id_1].size();
        }

        ssd += ssd_1 / den;
    }

    ssd_ag = ssd_total - ssd;

    // cout << "  ssd_ag: "<< ssd_ag << "\n";

    return ssd_ag;
}

double
LocusDivergence::haplotype_d_est(const Datum **d, const LocSum **s, vector<int> &pop_ids)
{
    //
    // Calculate D_est, fixation index, as described by
    //   Bird, et al., 2011, Detecting and measuring genetic differentiation
    //     +-Equation 11
    // and
    //   Jost, 2008, GST and its relatives do not measure differentiation, Molecular Ecology
    //     +- Equation 13, D_est_chao
    //
    map<string, double>            loc_haplotypes;
    map<int, map<string, double> > pop_haplotypes;
    map<int, double>               pop_totals;

    map<string, double>::iterator it;

    uint pop_cnt = pop_ids.size();

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int pop_id : pop_ids) {
        const Pop& pop = this->_mpopi->pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) {
                continue;
            } else if (d[i]->obshap.size() > 2) {
                continue;
            } else if (d[i]->obshap.size() == 1) {
                loc_haplotypes[d[i]->obshap[0]]         += 2;
                pop_haplotypes[pop_id][d[i]->obshap[0]] += 2;

            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    loc_haplotypes[d[i]->obshap[j]]++;
                    pop_haplotypes[pop_id][d[i]->obshap[j]]++;
                }
            }
        }

        for (it = pop_haplotypes[pop_id].begin(); it != pop_haplotypes[pop_id].end(); it++)
            pop_totals[pop_id] += it->second;
    }

    double x = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {

        double freq_sum_sq = 0.0;
        double freq_sq_sum = 0.0;
        for (int pop_id : pop_ids) {
            freq_sum_sq += (pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]);
            freq_sq_sum += pow((pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]), 2);
        }
        freq_sum_sq = pow(freq_sum_sq, 2);

        x += (freq_sum_sq - freq_sq_sum) / (pop_cnt - 1);
    }

    double y = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {
        for (int pop_id : pop_ids) {
            y += (pop_haplotypes[pop_id][it->first] * (pop_haplotypes[pop_id][it->first] - 1)) /
                (pop_totals[pop_id] * (pop_totals[pop_id] - 1));
        }
    }

    double d_est = 1.0 - (x / y);

    return d_est;
}

double
LocusDivergence::haplotype_dxy(const Datum **d, size_t loc_len, vector<int> &pop_ids)
{
    //
    // Calculate Dxy, fixation index, as described by
    //   Nei, 1987, Molecular Evolutionary Genetics, Chapter 10
    //     +-Equation 10.20 and Equation 5.3
    // and
    //   Cruickshank & Hahn, 2014, Reanalysis suggests that genomic islands of speciation are due
    //   to reduced diversity, not reduced gene flow. Molecular Ecology
    //     +- Box 1
    //
    map<string, double>            loc_haplotypes;
    map<int, map<string, double> > pop_haplotypes;
    map<int, double>               pop_totals;

    map<string, double>::iterator it;

    uint pop_cnt = pop_ids.size();

    assert(pop_cnt == 2);
    
    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int pop_id : pop_ids) {
        const Pop& pop = this->_mpopi->pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL)
                continue;

            assert(d[i]->obshap.size() == 2);

            for (uint j = 0; j < d[i]->obshap.size(); j++) {
                loc_haplotypes[d[i]->obshap[j]]++;
                pop_haplotypes[pop_id][d[i]->obshap[j]]++;
            }
        }

        for (it = pop_haplotypes[pop_id].begin(); it != pop_haplotypes[pop_id].end(); it++)
            pop_totals[pop_id] += it->second;
    }

    vector<string> haps;
    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++)
        haps.push_back(it->first);

    double popx_freq, popy_freq, nuc_diff;
    const char *p, *q;
    
    double dxy  = 0.0;
    uint   popx = pop_ids[0];
    uint   popy = pop_ids[1];
    for (uint i = 0; i < haps.size(); i++) {
        popx_freq = pop_haplotypes[popx][haps[i]] / pop_totals[popx];

        for (uint j = 0; j < haps.size(); j++) {
            popy_freq = pop_haplotypes[popy][haps[j]] / pop_totals[popy];
            
            nuc_diff  = 0;
            p = haps[i].c_str();
            q = haps[j].c_str();
            for (; *p != '\0'; p++, q++)
                nuc_diff += *p != *q ? 1 : 0;
            nuc_diff = nuc_diff / (double) loc_len;

            // Nei 1987, Equation 5.3
            nuc_diff = -1 * (3.0 / 4.0) * log(1 - ((4.0 / 3.0) * nuc_diff));

            // Nei 1987, Equation 10.20
            dxy += popx_freq * popy_freq * nuc_diff;
        }
    }

    return dxy;
}

bool
LocusDivergence::fixed_locus(const Datum **d, vector<int> &pop_ids)
{
    set<string>               loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;

    for (int pop_id : pop_ids) {
        const Pop& pop = this->_mpopi->pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) continue;

            if (d[i]->obshap.size() > 2) {
                continue;

            } else if (d[i]->obshap.size() == 1) {
                if (!uncalled_haplotype(d[i]->obshap[0])) {
                    loc_haplotypes.insert(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                }
            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    if (!uncalled_haplotype(d[i]->obshap[0])) {
                        loc_haplotypes.insert(d[i]->obshap[j]);
                        pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
                    }
                }
            }
        }
    }

    uint valid_pops = 0;

    for (int pop_id : pop_ids) {
        if (pop_haplotypes[pop_id].size() > 0)
            valid_pops++;
    }

    //
    // Check that more than one population has data for this locus.
    //
    if (valid_pops <= 1)
        return true;

    //
    // Check that there is more than one haplotype at this locus.
    //
    if (loc_haplotypes.size() == 1)
        return true;

    return false;
}

int
LocusDivergence::nuc_substitution_identity(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {

            if (haplotypes[i] == haplotypes[j])
                dist = 0.0;
            else
                dist = 1.0;

            hdists[i][j] = dist;
            hdists[j][i] = dist;
        }
    }

    return 0;
}

int
LocusDivergence::nuc_substitution_identity_max(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {
            hdists[i][j] = 1.0;
            hdists[j][i] = 1.0;
        }
    }

    return 0;
}

int
LocusDivergence::write_summary(string path)
{
    //
    // Write out the mean Fst measure of each pair of populations.
    //
    string file = path + ".fst_summary.tsv";
    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening fst summary output file '" << file << "'\n";
        exit(1);
    }

    //
    // Write out X-axis header.
    //
    for (auto& pop : this->_mpopi->pops())
        fh << "\t" << pop.name;
    fh << "\n";

    uint n = 0;
    for (uint i = 0; i < this->_mpopi->pops().size() - 1; i++) {
        fh << this->_mpopi->pops()[i].name;

        for (uint k = 0; k <= i; k++)
            fh << "\t";

        for (uint j = i + 1; j < this->_mpopi->pops().size(); j++) {
            fh << "\t" << this->_mean_fst[n] / this->_mean_fst_cnt[n];
            n++;
        }
        fh << "\n";
    }

    fh.close();

    file = path + ".phistats_summary.tsv";
    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening fst summary output file '" << file << "'\n";
        exit(1);
    }

    fh << "# Phi_st Means\n";
    //
    // Write out X-axis header.
    //
    for (auto& pop : this->_mpopi->pops())
        fh << "\t" << pop.name;
    fh << "\n";

    n = 0;
    for (uint i = 0; i < this->_mpopi->pops().size() - 1; i++) {
        fh << this->_mpopi->pops()[i].name;

        for (uint k = 0; k <= i; k++)
            fh << "\t";

        for (uint j = i + 1; j < this->_mpopi->pops().size(); j++) {
            fh << "\t" << this->_mean_phist[n] / this->_mean_phist_cnt[n];
            n++;
        }
        fh << "\n";
    }

    fh << "\n"
       << "# Fst' Means\n";
    //
    // Write out X-axis header.
    //
    for (auto& pop : this->_mpopi->pops())
        fh << "\t" << pop.name;
    fh << "\n";

    n = 0;
    for (uint i = 0; i < this->_mpopi->pops().size() - 1; i++) {
        fh << this->_mpopi->pops()[i].name;

        for (uint k = 0; k <= i; k++)
            fh << "\t";

        for (uint j = i + 1; j < this->_mpopi->pops().size(); j++) {
            fh << "\t" << this->_mean_fstp[n] / this->_mean_fstp_cnt[n];
            n++;
        }
        fh << "\n";
    }

    fh << "\n"
       << "# Dxy Means\n";
    //
    // Write out X-axis header.
    //
    for (auto& pop : this->_mpopi->pops())
        fh << "\t" << pop.name;
    fh << "\n";

    n = 0;
    for (uint i = 0; i < this->_mpopi->pops().size() - 1; i++) {
        fh << this->_mpopi->pops()[i].name;

        for (uint k = 0; k <= i; k++)
            fh << "\t";

        for (uint j = i + 1; j < this->_mpopi->pops().size(); j++) {
            fh << "\t" << this->_mean_dxy[n] / this->_mean_dxy_cnt[n];
            n++;
        }
        fh << "\n";
    }

    fh.close();

    cout << "\nPopulation pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):\n";
    n = 0;
    for (uint i = 0; i < this->_mpopi->pops().size() - 1; i++) {

        for (uint j = i + 1; j < this->_mpopi->pops().size(); j++) {
            cout << "  " << this->_mpopi->pops()[i].name << "-"
                 << this->_mpopi->pops()[j].name
                 << ": mean Fst: "    << this->_mean_fst[n]   / this->_mean_fst_cnt[n]
                 << "; mean Phi_st: " << this->_mean_phist[n] / this->_mean_phist_cnt[n]
                 << "; mean Fst': "   << this->_mean_fstp[n]  / this->_mean_fstp_cnt[n]
                 << "; mean Dxy: "    << this->_mean_dxy[n]   / this->_mean_dxy_cnt[n]
                 << "\n";
            n++;
        }
    }
    cout << "\n";

    return 0;
}
