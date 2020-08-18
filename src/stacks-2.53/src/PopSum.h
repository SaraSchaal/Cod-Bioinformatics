// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2019, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __POPSUM_H__
#define __POPSUM_H__

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include <cstdint>
#include <cmath>

#include "stacks.h"
#include "locus.h"
#include "PopMap.h"
#include "MetaPopInfo.h"
#include "Hwp.h"

extern bool      calc_hwp;
extern bool      log_fst_comp;
extern double    minor_allele_freq;
extern corr_type fst_correction;
extern double    p_value_cutoff;

const uint PopStatSize = 6;

class PopStat {
public:
    int      loc_id;
    int      bp;
    bool     fixed;
    double   alleles;    // Number of alleles sampled at this location.
    uint16_t snp_cnt;    // Number of SNPs in kernel-smoothed window centered on this SNP.
    double   stat[PopStatSize];
    mutable double smoothed[PopStatSize];
    double   bs[PopStatSize];

    PopStat() {
        this->loc_id   = 0;
        this->bp       = 0;
        this->fixed    = false;
        this->alleles  = 0.0;
        this->snp_cnt  = 0;

        for (uint i = 0; i < PopStatSize; i++) {
            this->stat[i]     = 0.0;
            this->smoothed[i] = 0.0;
            this->bs[i]       = 0.0;
        }
    }
    int reset() {
        this->loc_id   = 0;
        this->bp       = 0;
        this->fixed    = false;
        this->alleles  = 0.0;
        this->snp_cnt  = 0;

        for (uint i = 0; i < PopStatSize; i++) {
            this->stat[i]     = 0.0;
            this->smoothed[i] = 0.0;
            this->bs[i]       = 0.0;
        }
        return 0;
    }
    virtual ~PopStat() {}
};

class HapStat: public PopStat {
    // PopStat[0]: Phi_st
    // PopStat[1]: Phi_ct
    // PopStat[2]: Phi_sc
    // PopStat[3]: Fst'
    // PopStat[4]: D_est
    // PopStat[5]: Dxy
public:
    uint    pop_1;
    uint    pop_2;
    double *comp;
    uint    popcnt;

    HapStat(): PopStat() {
        comp = NULL;
        popcnt = uint(-1);
    }
    ~HapStat() {
        if (this->comp != NULL)
            delete [] comp;
    }
};

class LocStat: public PopStat {
    // PopStat[0]: gene diversity
    // PopStat[1]: haplotype diversity (Pi)
    // PopStat[2]: Hardy-Weinberg proportions p-value
    // PopStat[3]: Hardy-Weinberg proportions p-value standard error.
public:
    uint     hap_cnt; // Number of unique haplotypes at this locus.
    string   hap_str; // Human-readable string of haplotype counts.

    LocStat(): PopStat() {
        this->hap_cnt = 0;
    }
    ~LocStat() {};
};

class PopPair: public PopStat {
    // PopStat[0]: corrected Fst, (by p-value or Bonferroni p-value).
    // PopStat[1]: corrected AMOVA Fst
public:
    uint    pop_1;
    uint    pop_2;
    uint    col;
    double  pi;
    double  fst;
    double  fet_p;      // Fisher's Exact Test p-value.
    double  fet_or;     // Fisher's exact test odds ratio.
    double  or_se;      // Fisher's exact test odds ratio standard error.
    double  lod;        // base 10 logarithm of odds score.
    double  ci_low;     // Fisher's exact test lower confidence interval.
    double  ci_high;    // Fisher's exact test higher confidence interval.
    double  amova_fst;  // AMOVA Fst method, from Weir, Genetic Data Analysis II .
    double *comp;

    PopPair() {
        pop_1     = 0;
        pop_2     = 0;
        col       = 0;
        pi        = 0.0;
        fst       = 0.0;
        fet_p     = 0.0;
        fet_or    = 0.0;
        or_se     = 0.0;
        lod       = 0.0;
        ci_low    = 0.0;
        ci_high   = 0.0;
        amova_fst = 0.0;
        comp      = NULL;
    }
    ~PopPair() {
        if (this->comp != NULL)
            delete [] comp;
    }
};

class SumStat: public PopStat {
    // PopStat[0]: pi
    // PopStat[1]: fis
    // PopStat[2]: HWE p-value
public:
    bool    incompatible_site; // Deprecated? c.f. Note from 2018-12-17 in LocPopSum::sum_pops().
    bool    filtered_site;
    double  num_indv;
    char    p_nuc;
    char    q_nuc;
    double  p;
    double  obs_het;
    double  obs_hom;
    double  exp_het;
    double  exp_hom;
    double &pi;

    SumStat(): PopStat(), pi(this->stat[0]) {
        num_indv  = 0.0;
        p         = 0.0;
        p_nuc     = 0;
        q_nuc     = 0;
        obs_het   = 0.0;
        obs_hom   = 0.0;
        exp_het   = 0.0;
        exp_hom   = 0.0;
        snp_cnt   = 0;
        incompatible_site = false;
        filtered_site     = false;

        this->stat[2] = 1.0; // The default value for the HWE p-value should be 1.0, not 0.
    }
    int reset() {
        PopStat::reset();
        num_indv  = 0.0;
        p         = 0.0;
        p_nuc     = 0;
        q_nuc     = 0;
        obs_het   = 0.0;
        obs_hom   = 0.0;
        exp_het   = 0.0;
        exp_hom   = 0.0;
        snp_cnt   = 0;
        incompatible_site = false;
        filtered_site     = false;
        return 0;
    }
    ~SumStat() {assert(incompatible_site == false);} // c.f. Note from 2018-12-17 in LocPopSum::sum_pops().
};

class LocSum {
public:
    SumStat *nucs;    // Array containing summary statistics for
                      // each nucleotide position at this locus.
    LocSum(int len) {
        this->nucs = new SumStat[len];
    }
    ~LocSum() {
        delete [] this->nucs;
    }
};

class NucTally {
public:
    int      loc_id;
    int      bp;
    uint16_t col;
    uint16_t num_indv;
    uint16_t pop_cnt;
    uint16_t allele_cnt;
    char     p_allele;
    char     q_allele;
    double   p_freq;
    double   obs_het;
    bool     fixed;
    int      priv_allele;

    NucTally() {
        loc_id      = 0;
        bp          = 0;
        col         = 0;
        num_indv    = 0;
        pop_cnt     = 0;
        allele_cnt  = 0;
        p_allele    = 0;
        q_allele    = 0;
        p_freq      = 0.0;
        obs_het     = 0.0;
        priv_allele = -1;
        fixed       = true;
    }
    int reset() {
        loc_id      = 0;
        bp          = 0;
        col         = 0;
        num_indv    = 0;
        pop_cnt     = 0;
        allele_cnt  = 0;
        p_allele    = 0;
        q_allele    = 0;
        p_freq      = 0.0;
        obs_het     = 0.0;
        priv_allele = -1;
        fixed       = true;
        return 0;
    }
};

class LocTally {
public:
    NucTally *nucs;

    LocTally(int len)  {
        this->nucs = new NucTally[len];
    }
    ~LocTally() {
        delete [] this->nucs;
    }
};

//
// per-Locus Population-level summary
//
class LocPopSum {
    size_t    _pop_cnt;
    LocSum  **_per_pop;
    LocTally *_meta_pop;
    LocStat **_hapstats_per_pop;

public:
    LocPopSum(size_t cloc_len, const MetaPopInfo& mpopi);
    ~LocPopSum();

    int             sum_pops(const CSLocus *, Datum const*const*, const MetaPopInfo&, bool, ostream &);
    int             tally_metapop(const CSLocus *);
    int             calc_hapstats(const CSLocus *, const Datum **, const MetaPopInfo&);
    const LocSum  **all_pops()                         { return (const LocSum **) this->_per_pop; }
    const LocSum   *per_pop(size_t pop_index)          const { return this->_per_pop[pop_index]; }
    const LocTally *meta_pop()                         const { return this->_meta_pop; }
    const LocStat  *hapstats_per_pop(size_t pop_index) const { return this->_hapstats_per_pop[pop_index]; }
    size_t          pop_cnt()                          const { return this->_pop_cnt; }
    static double   pi(double, double, double);
    static double   binomial_coeff(double, double);
    double          hwe(double, double, double, double, double, double);

private:
    int      tally_heterozygous_pos(const CSLocus *, Datum const*const*, LocSum *, int, int, uint, uint);
    int      tally_fixed_pos(const CSLocus *, Datum const*const*, LocSum *, int, uint, uint);
    int      tally_ref_alleles(int, uint16_t &, char &, char &, uint16_t &, uint16_t &);
    int      tally_observed_haplotypes(const vector<char *> &, int);
    LocStat *haplotype_diversity(int, int, Datum const*const*);
    double   log_hwp_pr(double, double, double, double, double, double);
};

struct LocBin {
    size_t     sample_cnt;
    CSLocus   *cloc;
    Datum    **d;
    LocPopSum *s;

    LocBin(size_t cnt): sample_cnt(cnt), cloc(NULL), d(NULL), s(NULL) {}
    ~LocBin() {
        if (this->cloc != NULL) delete cloc;
        if (this->d    != NULL) {
            for (uint i = 0; i < sample_cnt; i++)
                if (d[i] != NULL) delete d[i];
            delete [] d;
        }
        if (this->s    != NULL) delete s;
    }
};

//
// per-Locus class for calculating divergence values such as Fst.
//
class LocusDivergence {
    const MetaPopInfo         *_mpopi;
    vector<vector<PopPair **>> _snps;
    vector<vector<HapStat *>>  _haplotypes;
    vector<HapStat *>          _metapop_haplotypes;

    vector<double> _mean_fst;
    vector<double> _mean_fst_cnt;
    vector<double> _mean_phist;
    vector<double> _mean_phist_cnt;
    vector<double> _mean_fstp;
    vector<double> _mean_fstp_cnt;
    vector<double> _mean_dxy;
    vector<double> _mean_dxy_cnt;

public:
    LocusDivergence(const MetaPopInfo *mpopi);
    int clear(const vector<LocBin *> &loci);
    int snp_divergence(const vector<LocBin *> &loci);
    int haplotype_divergence_pairwise(const vector<LocBin *> &loci);
    int haplotype_divergence(const vector<LocBin *> &loci);

    int write_summary(string);

    vector<vector<PopPair **>>& snp_values()       { return this->_snps; }
    vector<vector<HapStat *>>&  haplotype_values() { return this->_haplotypes; }
    vector<HapStat *>&  metapop_haplotype_values() { return this->_metapop_haplotypes; }

private:
    //
    // SNP-level F statistics.
    //
    PopPair *Fst(const CSLocus *, const LocPopSum *, int, int, int);
    int      fishers_exact_test(PopPair *, double, double, double, double);

    //
    // Haplotype-level F statistics
    //
    double   haplotype_d_est(const Datum **, const LocSum **, vector<int> &);
    HapStat *haplotype_amova(const Datum **, const LocSum **, vector<int> &);
    double   amova_ssd_total(vector<string> &, map<string, int> &, double **);
    double   amova_ssd_wp(vector<int> &, map<int, vector<int>> &, map<string, int> &, map<int, vector<string>> &, double **);
    double   amova_ssd_ap_wg(vector<int> &, map<int, vector<int>> &, map<string, int> &, map<int, vector<string>> &, double **, double **);
    double   amova_ssd_ag(vector<int> &, map<int, vector<int>> &, map<string, int> &, map<int, vector<string>> &, double **, double);

    //
    // Haplotype-level Dxy
    //
    double   haplotype_dxy(const Datum **, size_t,  vector<int> &);

    bool     fixed_locus(const Datum **, vector<int> &);
    int      nuc_substitution_identity(map<string, int> &, double **);
    int      nuc_substitution_identity_max(map<string, int> &, double **);
};

//
// Utility functions for summary/haplotype statistics.
//
bool     uncalled_haplotype(const char *);
double   count_haplotypes_at_locus(int, int, Datum const*const*, map<string, double> &);
int      nuc_substitution_dist(map<string, int> &, double **);


#endif // __POPSUM_H__
