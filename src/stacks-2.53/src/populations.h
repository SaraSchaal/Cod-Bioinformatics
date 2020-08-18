// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2019, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __POPULATIONS_H__
#define __POPULATIONS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cctype>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <set>

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "renz.h"
#include "PopMap.h"
#include "PopSum.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "mapping_utils.h"
#include "sql_utilities.h"
#include "genotype_dictionaries.h"
#include "ordered.h"
#include "smoothing.h"
#include "bootstrap.h"
#include "MetaPopInfo.h"
#include "gzFasta.h"
#include "locus_readers.h"

enum bs_type   {bs_exact, bs_approx, bs_none};
enum merget    {merge_sink, merge_src};
enum phaset    {merge_failure, simple_merge, complex_phase, nomapping_fail, multimapping_fail, multiple_fails};
enum class InputMode {stacks, stacks2, vcf};

//
// Class for accumulating summary statistics as we read each batch of data.
// After all loci are processed, output the summary statistics summary.
//
class SumStatsSummary {
    size_t  _pop_cnt;
    size_t *_private_cnt;
    size_t *_sig_hwe_dev;
    double *_num_indv_mean,         *_p_mean,         *_obs_het_mean,         *_obs_hom_mean,         *_exp_het_mean,         *_exp_hom_mean,         *_pi_mean,         *_fis_mean;
    double *_num_indv_acc_mean,     *_p_acc_mean,     *_obs_het_acc_mean,     *_obs_hom_acc_mean,     *_exp_het_acc_mean,     *_exp_hom_acc_mean,     *_pi_acc_mean,     *_fis_acc_mean;
    double *_num_indv_var,          *_p_var,          *_obs_het_var,          *_obs_hom_var,          *_exp_het_var,          *_exp_hom_var,          *_pi_var,          *_fis_var;
    double *_num_indv_mean_all,     *_p_mean_all,     *_obs_het_mean_all,     *_obs_hom_mean_all,     *_exp_het_mean_all,     *_exp_hom_mean_all,     *_pi_mean_all,     *_fis_mean_all;
    double *_num_indv_acc_mean_all, *_p_acc_mean_all, *_obs_het_acc_mean_all, *_obs_hom_acc_mean_all, *_exp_het_acc_mean_all, *_exp_hom_acc_mean_all, *_pi_acc_mean_all, *_fis_acc_mean_all;
    double *_num_indv_var_all,      *_p_var_all,      *_obs_het_var_all,      *_obs_hom_var_all,      *_exp_het_var_all,      *_exp_hom_var_all,      *_pi_var_all,      *_fis_var_all;
    double *_n, *_n_all, *_var_sites;
    double *_sq_n, *_sq_n_all;
    double _locus_len_mean, _locus_len_acc_mean, _locus_len_var, _locus_len_mean_all, _locus_len_acc_mean_all, _locus_len_var_all;
    double _locus_gt_sites_mean, _locus_gt_sites_acc_mean, _locus_gt_sites_var;
    double _overlap_mean, _overlap_acc_mean, _overlap_var;
    double _locus_n, _locus_pe_ctg_n, _locus_overlap_n;

public:
    SumStatsSummary(size_t pop_cnt);
    ~SumStatsSummary();
    int accumulate(const vector<LocBin *> &);
    int final_calculation();
    int write_results();

private:
    double online_variance(double x, double &acc_mean, double n);
};

//
// Class for storing distributions of data related to catalog loci.
//
class CatalogDists {
    map<size_t, size_t> _pre_valid,  _pre_absent,  _pre_confounded;  // Prior to any filtering.
    map<size_t, size_t> _post_valid, _post_absent, _post_confounded; // After filtering.
    map<size_t, size_t> _pre_snps_per_loc, _post_snps_per_loc;

public:
    int accumulate_pre_filtering(const size_t, const CSLocus *);
    int accumulate(const vector<LocBin *> &);
    int write_results(ostream &log_fh);
};

//
// Class for smoothing various calculated population statistics for each batch of loci.
//
class LocusSmoothing {
    const MetaPopInfo    *_mpopi; // Population Map
    KSmooth<SumStat>     *_ks_ss;
    OSumStat<SumStat>    *_ord_ss;
    KSmooth<LocStat>     *_ks_ls;
    OHaplotypes<LocStat> *_ord_ls;
    KSmooth<HapStat>     *_ks_hs;
    OHaplotypes<HapStat> *_ord_hs;
    KSmooth<PopPair>     *_ks_pp;
    OPopPair<PopPair>    *_ord_pp;

public:
    LocusSmoothing(const MetaPopInfo *mpopi, ofstream &log_fh) {
        this->_mpopi  = mpopi;
        this->_ks_ss  = new KSmooth<SumStat>(2);
        this->_ord_ss = new OSumStat<SumStat>(log_fh);
        this->_ks_ls  = new KSmooth<LocStat>(2);
        this->_ord_ls = new OHaplotypes<LocStat>();
        this->_ks_hs  = new KSmooth<HapStat>(6);
        this->_ord_hs = new OHaplotypes<HapStat>();
        this->_ks_pp  = new KSmooth<PopPair>(2);
        this->_ord_pp = new OPopPair<PopPair>(log_fh);
    }
    ~LocusSmoothing() {
        delete this->_ks_ss;
        delete this->_ord_ss;
        delete this->_ks_ls;
        delete this->_ord_ls;
        delete this->_ks_hs;
        delete this->_ord_hs;
        delete this->_ks_pp;
        delete this->_ord_pp;
    }

    int    snpstats(const vector<LocBin *> &, ofstream &);
    int    hapstats(const vector<LocBin *> &, ofstream &);
    int    snp_divergence(const vector<LocBin *> &, const vector<vector<PopPair **>> &, ofstream &);
    int    hap_divergence(const vector<LocBin *> &, const vector<vector<HapStat *>>  &, const vector<HapStat *> &, ofstream &);
};

//
// Class for filtering whole loci based on the sample and population limits (-r, -p).
//
class LocusFilter {
public:
    LocusFilter() {
        this->_pop_cnt    = 0;
        this->_sample_cnt = 0;
        this->_pop_order  = NULL; // The array order of each population.
        this->_samples    = NULL; // Which population each sample belongs to.
        this->_pop_tot    = NULL; // The total number of samples in each population.
        this->_filtered_loci       = 0;
        this->_total_loci          = 0;
        this->_seen_loci           = 0;
        this->_batch_filtered_loci = 0;
        this->_batch_total_loci    = 0;
        this->_batch_seen_loci     = 0;
        this->_filtered_sites      = 0;
        this->_total_sites         = 0;
        this->_variant_sites       = 0;
    }
    LocusFilter(MetaPopInfo *mpopi) {
        assert(mpopi != NULL);
        this->init(mpopi);
    }
    ~LocusFilter() {
        if (this->_pop_order != NULL)
            delete [] this->_pop_order;
        if (this->_samples != NULL)
            delete [] this->_samples;
        if (this->_pop_tot != NULL)
            delete [] this->_pop_tot;
    }
    // Remove the copy and move constructors.
    LocusFilter(const LocusFilter &) = delete;
    LocusFilter(LocusFilter &&) = delete;

    int    load_whitelist(string);
    int    load_blacklist(string);

    void   init(MetaPopInfo *mpopi);
    bool   whitelist_filter(size_t locus_id);
    bool   blacklist_filter(size_t locus_id);
    void   whitelist_snp_filter(LocBin& loc) const;
    bool   apply_filters_stacks(LocBin& loc, ostream& log_fh, const MetaPopInfo& mpopi);
    bool   apply_filters_external(LocBin& loc, ostream& log_fh, const MetaPopInfo& mpopi);
    bool   filter(const MetaPopInfo *mpopi, Datum **d, CSLocus *cloc);
    void   filter_snps(LocBin& loc, const MetaPopInfo& mpopi, ostream &log_fh);
    void   filter_haps(LocBin& loc, const MetaPopInfo& mpopi, ostream &log_fh);
    void   gt_depth_filter(Datum **d, const CSLocus *cloc);
    void   keep_single_snp(CSLocus* cloc, Datum** d, size_t n_samples, const LocTally* t) const;
    void   keep_random_snp(CSLocus* cloc, Datum** d, size_t n_samples, const LocTally* t) const;

    size_t filtered()       const { return this->_filtered_loci; }
    size_t total()          const { return this->_total_loci; }
    size_t seen()           const { return this->_seen_loci; }
    size_t batch_filtered() const { return this->_batch_filtered_loci; }
    size_t batch_seen()     const { return this->_batch_seen_loci; }
    size_t batch_total()    const { return this->_batch_total_loci; }
    size_t filtered_sites() const { return this->_filtered_sites; }
    size_t total_sites()    const { return this->_total_sites; }
    size_t variant_sites()  const { return this->_variant_sites; }
    void   locus_seen();
    void   locus_unsee();
    void   keep_locus(LocBin *);
    void   batch_clear();

    const set<int>&            blacklist() { return this->_blacklist; }
    const map<int, set<int>>&  whitelist() { return this->_whitelist; }

private:
    static void erase_snp(CSLocus *cloc, Datum **d, size_t n_samples, size_t snp_index);

    size_t  _pop_cnt;
    size_t  _sample_cnt;
    size_t *_pop_order;
    size_t *_samples;
    size_t *_pop_tot;
    size_t  _filtered_loci;
    size_t  _total_loci;
    size_t  _seen_loci;
    size_t  _batch_filtered_loci;
    size_t  _batch_total_loci;
    size_t  _batch_seen_loci;
    size_t  _filtered_sites;
    size_t  _total_sites;
    size_t  _variant_sites;

    set<int>           _blacklist;
    map<int, set<int>> _whitelist;
};

//
// BatchLocusProcessor
// ----------
// Class for processing loci in batches, or per chromosome.
//
class BatchLocusProcessor {
public:
    BatchLocusProcessor():
        _input_mode(InputMode::stacks2), _user_supplied_whitelist(false), _batch_size(0), _batch_num(0),
        _mpopi(NULL), _next_loc(NULL), _unordered_bp(0) {}
    BatchLocusProcessor(InputMode mode, size_t batch_size, MetaPopInfo *popi):
        _input_mode(mode), _user_supplied_whitelist(false), _batch_size(batch_size), _batch_num(0),
        _mpopi(popi), _next_loc(NULL), _unordered_bp(0) {}
    BatchLocusProcessor(InputMode mode, size_t batch_size):
        _input_mode(mode), _user_supplied_whitelist(false), _batch_size(batch_size), _batch_num(0),
        _mpopi(NULL), _next_loc(NULL), _unordered_bp(0) {}
    ~BatchLocusProcessor() {
        for (uint i = 0; i < this->_loci.size(); i++)
            delete this->_loci[i];
        delete [] this->_sig_hwe_dev;
    };

    int            init(string, string);
    size_t         next_batch(ostream &);
    size_t         next_batch_number() { return this->_batch_num + 1; }
    int            summarize(ostream &);
    int            hapstats(ostream &);
    int            write_distributions(ostream &log_fh) { return this->_dists.write_results(log_fh); }

    MetaPopInfo*   pop_info()      { return this->_mpopi; }
    int            pop_info(MetaPopInfo *popi) { this->_mpopi = popi; return 0; }
    VcfParser&     vcf_reader()    { return this->_vcf_parser; }
    GzFasta&       fasta_reader()  { return this->_fasta_reader; }
    VcfCLocReader& cloc_reader()   { return this->_cloc_reader; }
    size_t         batch_size()    { return this->_batch_size; }
    size_t         batch_size(size_t bsize) { this->_batch_size = bsize; return bsize; }
    void           report_locus_overlap(size_t&, size_t&, ofstream* = NULL);

    const LocusFilter&      filter() { return this->_loc_filter; }
    const vector<LocBin *>& loci()   { return this->_loci; }
    const string&           chr()    { return this->_chr; }
    const CatalogDists&     dists()  { return this->_dists; }

    // Per-population haplotype counters
    size_t           *_sig_hwe_dev;

private:
    InputMode    _input_mode;
    bool         _user_supplied_whitelist;
    size_t       _batch_size; // Number of loci to process at a time.
    size_t       _batch_num;  // Counter for how many batches we have processed.
    MetaPopInfo *_mpopi;      // Population Map

    // Parsers
    VcfParser     _vcf_parser;
    VcfCLocReader _cloc_reader;
    GzFasta       _fasta_reader;
    vector<size_t> _samples_vcf_to_mpopi;

    // Data stores
    vector<LocBin *>  _loci;
    LocBin           *_next_loc;
    string            _chr;
    CatalogDists      _dists;

    // Counters for external VCF
    size_t            _total_ext_vcf;
    vector<size_t>    _skipped_notsnp;
    vector<size_t>    _skipped_notbinarysnp;
    vector<size_t>    _skipped_filter;

    // Controls for which loci are loaded
    LocusFilter _loc_filter;

    size_t _unordered_bp;

private:
    int    init_external_loci(string, string);
    int    init_stacks_loci(string, string);
    void   batch_clear();
    size_t next_batch_external_loci(ostream &);
    size_t next_batch_stacks_loci(ostream &);
};

void    help( void );
void    version( void );
int     parse_command_line(int, char**);
void    output_parameters(ostream &);
void    open_log(ofstream &);
int     build_file_list();
bool    check_population_map_for_cross(const MetaPopInfo *);
int     load_marker_list(string, set<int> &);
int     load_marker_column_list(string, map<int, set<int> > &);
//int     merge_shared_cutsite_loci(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &, ofstream &);
phaset  merge_and_phase_loci(PopMap<CSLocus> *, CSLocus *, CSLocus *, set<int> &, ofstream &);
int     merge_datums(int, int, Datum **, Datum **, set<string> &, int);
int     merge_csloci(CSLocus *, CSLocus *, set<string> &);
int     tabulate_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *);
int     tabulate_locus_haplotypes(CSLocus *, Datum **, int);
int     create_genotype_map(CSLocus *, Datum **, int);
int     call_population_genotypes(CSLocus *, Datum **, int);
int     translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &); // This function doesn't exist (March 24, 2016)
int     correct_fst_bonferroni_win(vector<PopPair *> &);
int     bootstrap_fst_approximate_dist(vector<double> &, vector<int>  &, double *, int *, map<int, vector<double> > &); // not used (March 23, 2016)
int     bootstrap_popstats_approximate_dist(vector<double> &, vector<double> &, vector<int>  &, double *, int *, int, map<int, vector<double> > &, map<int, vector<double> > &); // not used (March 23, 2016)
double  bootstrap_approximate_pval(int, double, map<int, vector<double> > &);

bool hap_compare(const pair<string,int>&, const pair<string,int>&);

void vcfcomp_simplify_pmap (map<int, CSLocus*>& catalog, PopMap<CSLocus>* pmap);

#endif // __POPULATIONS_H__
