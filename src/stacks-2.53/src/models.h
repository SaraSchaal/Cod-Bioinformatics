// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2012, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __MODELS_H__
#define __MODELS_H__

#include "constants.h"
#include "utils.h"
#include "DNASeq4.h"
#include "stacks.h"
#include "locus.h"
#include "mstack.h"

//
// Possible models for calling nucleotide positions as fixed or variable
//
enum modelt {fixed, snp, bounded, marukihigh, marukilow};

//
// For use with the multinomial model to call fixed nucleotides.
//
extern int barcode_size;
extern double barcode_err_freq;

extern double heterozygote_limit;
extern double homozygote_limit;
extern double bound_low;  // For the bounded-snp model.
extern double bound_high; // For the bounded-snp model.
extern double p_freq;     // For the fixed model.

extern const std::array<double,40> chisq1ddf_gq;


bool lrtest(double lnl_althyp, double lnl_nullhyp, double threshold); // Where threshold is a value in the Chi2 distribution.
long vcf_lnl2gq(double lnl_gt, double lnl_gt_alt); // Returns VCF's GQ field `floor(-10 * log10(p-value))`

double qchisq(double alpha, size_t df);
modelt parse_model_type(const string& arg);
void set_model_thresholds(double alpha);
void report_model(ostream& os, modelt model_type);
void report_alpha(ostream& os, double alpha);
string to_string(modelt model_type);

snp_type call_snp (double l_ratio);
snp_type call_snp (double lnl_hom, double lnl_het);

double   lnl_multinomial_model_hom (double total, double n1);
double   lnl_multinomial_model_het (double total, double n1n2);
double   lr_multinomial_model         (double nuc_1, double nuc_2, double nuc_3, double nuc_4);
double   lr_bounded_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4);

void call_bounded_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
void call_bounded_multinomial_snp(Locus *, int, map<char, int> &);
void call_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
void call_multinomial_snp(Locus *, int, map<char, int> &);
void call_multinomial_fixed(MergedStack *, int, map<char, int> &);

double   heterozygous_likelihood(int, map<char, int> &);
double   homozygous_likelihood(int, map<char, int> &);

class SampleCall {
    GtLiks lnls_;
    // The genotype call and the corresponding nucleotides.
    // hom {nt, Nt2()} | het {min_nt, max_nt} | unk {Nt2(), Nt2()}
    // For hets, the two nucleotides are sorted lexically (A<C<G<T).
    snp_type call_;
    array<Nt2,2> nts_;
    long gq_;
public:
    SampleCall() : lnls_(), call_(snp_type_unk), nts_{{Nt2(),Nt2()}}, gq_(-1) {}

    const GtLiks& lnls() const {return lnls_;}
          GtLiks& lnls()       {return lnls_;}

    snp_type call() const {return call_;}
    Nt2 nt0() const {assert(call_==snp_type_hom || call_==snp_type_het); return nts_[0];}
    Nt2 nt1() const {assert(call_==snp_type_het); return nts_[1];}
    array<Nt2,2> nts() const {assert(call_==snp_type_het); return nts_;}
    long gq() const {assert(call_==snp_type_hom || call_==snp_type_het); return gq_;}

    void set_call(snp_type c, Nt2 rank0_nt, Nt2 rank1_nt, long gq);
    void discard() {call_ = snp_type_discarded; nts_= {{Nt2(), Nt2()}}; gq_ = -1;}

    // For debugging.
    friend ostream& operator<<(ostream& os, const SampleCall& sc);
};

class SiteCall {
    map<Nt2,double> alleles_;
    long snp_qual_;
    vector<SampleCall> sample_calls_; // Empty if alleles_.size() < 2.
public:
    SiteCall(
            map<Nt2,double>&& alleles,
            long snp_qual,
            vector<SampleCall>&& sample_calls
        )
            : alleles_(move(alleles))
            , snp_qual_(snp_qual)
            , sample_calls_(move(sample_calls))
        {}
    SiteCall()
        : SiteCall({}, -1, {}) {}
    SiteCall(Nt2 fixed_nt)
        : SiteCall({{fixed_nt, 1.0}}, -1, {}) {}
    SiteCall(map<Nt2,double>&& alleles, vector<SampleCall>&& sample_calls)
        : SiteCall(move(alleles), -1, move(sample_calls)) {}

    const map<Nt2,double>& alleles() const {return alleles_;}
    long snp_qual() const {return snp_qual_;}
    const vector<SampleCall>& sample_calls() const {return sample_calls_;}

    Nt2 most_frequent_allele() const;

    void filter_mac(size_t min_mac);
    void discard_sample(size_t sample_i)
        {assert(!sample_calls_.empty()); sample_calls_[sample_i].discard();}

    static Counts<Nt2> tally_allele_counts(const vector<SampleCall>& spldata);
    static map<Nt2,double> tally_allele_freqs(const vector<SampleCall>& spldata);

    // For debugging.
    void print(ostream& os, const SiteCounts& depths);
};

class Model {
public:
    virtual ~Model() {}
    virtual SiteCall call(const SiteCounts& depths) const = 0;
    virtual void print(ostream& os) const = 0;
    friend ostream& operator<< (ostream& os, const Model& m) {m.print(os); return os;}
};

//
// MultinomialModel: the standard Stacks v.1 model described in Hohenloe2010.
//
class MultinomialModel : public Model {
    double alpha_;
public:
    MultinomialModel(double gt_alpha) : alpha_(gt_alpha) {set_model_thresholds(alpha_);}
    SiteCall call(const SiteCounts& depths) const;
    void print(ostream& os) const
        {os << to_string(modelt::snp) << " (alpha: "  << alpha_ << ")";}
};

//
// MarukiHighModel: the model of Maruki & Lynch (2017) for high-coverage data.
//
class MarukiHighModel : public Model {
    double gt_alpha_;
    double gt_threshold_;
    double var_alpha_;
    double var_threshold_;
    double calc_hom_lnl(double n, double n1) const;
    double calc_het_lnl(double n, double n1n2) const;
public:
    MarukiHighModel(double gt_alpha, double var_alpha)
        : gt_alpha_(gt_alpha), gt_threshold_(qchisq(gt_alpha_,1)),
          var_alpha_(var_alpha), var_threshold_(qchisq(var_alpha_,1))
        {}
    SiteCall call(const SiteCounts& depths) const;
    void print(ostream& os) const
        {os << to_string(modelt::marukihigh) << " (var_alpha: "  << var_alpha_ << ", gt_alpha: " << gt_alpha_ << ")";}
};

//
// MarukiLowModel: the model of Maruki & Lynch (2015,2017) for low-coverage data.
//
class MarukiLowModel : public Model {
    struct LikData {
        bool has_data;
        double lnl_MM;
        double lnl_Mm;
        double lnl_mm;
        double l_MM;
        double l_Mm;
        double l_mm;
        LikData() : has_data(false), lnl_MM(0.0), lnl_Mm(0.0), lnl_mm(0.0), l_MM(1.0), l_Mm(1.0), l_mm(1.0) {}
        LikData(double lnl_MM_, double lnl_Mm_, double lnl_mm_)
            : has_data(true),
              lnl_MM(lnl_MM_), lnl_Mm(lnl_Mm_), lnl_mm(lnl_mm_),
              l_MM(exp(lnl_MM)), l_Mm(exp(lnl_Mm)), l_mm(exp(lnl_mm))
            {assert(std::isfinite(lnl_MM) && std::isfinite(lnl_Mm) && std::isfinite(lnl_mm));}
    };

    double gt_alpha_;
    double gt_threshold_;
    double var_alpha_;
    double var_threshold_;
    double calc_fixed_lnl(double n_tot, double n_M_tot) const;
    double calc_dimorph_lnl(double freq_MM, double freq_Mm, double freq_mm, const vector<LikData>& liks) const;
    double calc_ln_weighted_sum(double freq_MM, double freq_Mm, double freq_mm, const LikData& s_liks) const;
    double calc_ln_weighted_sum_safe(double freq_MM, double freq_Mm, double freq_mm, const LikData& s_liks) const;

    mutable size_t n_wsum_tot_;
    mutable size_t n_wsum_underflows_;
    mutable size_t n_called_sites_;
    mutable double sum_site_err_rates_;

public:
    MarukiLowModel(double gt_alpha, double var_alpha)
        : gt_alpha_(gt_alpha), gt_threshold_(qchisq(gt_alpha_,1)),
          var_alpha_(var_alpha), var_threshold_(qchisq(var_alpha_,2)), // df=2
          n_wsum_tot_(0), n_wsum_underflows_(0), n_called_sites_(0), sum_site_err_rates_(0.0)
        {}

    SiteCall call(const SiteCounts& depths) const;
    void print(ostream& os) const
        {os << to_string(modelt::marukilow) << " (var_alpha: "  << var_alpha_ << ", gt_alpha: " << gt_alpha_ << ")";}
    size_t n_wsum_tot() const {return n_wsum_tot_;}
    size_t n_wsum_underflows() const {return n_wsum_underflows_;}
    double mean_err_rate() const {return sum_site_err_rates_ / n_called_sites_;}
};

//
// ==================
// Inline definitions
// ==================
//

inline
bool lrtest(double lnl_althyp, double lnl_nullhyp, double threshold) {
    return 2.0 * (lnl_althyp - lnl_nullhyp) > threshold;
}

inline
long vcf_lnl2gq(double lnl_gt, double lnl_gt_alt) {
    double lr = 2.0 * (lnl_gt - lnl_gt_alt);
    auto itr = std::upper_bound(chisq1ddf_gq.begin(), chisq1ddf_gq.end(), lr);
    return itr - chisq1ddf_gq.begin();
}

inline
snp_type call_snp (double l_ratio) {
    if (l_ratio <= heterozygote_limit)
        return snp_type_het;
    else if (l_ratio >= homozygote_limit)
        return snp_type_hom;
    else
        return snp_type_unk;
}

inline
snp_type call_snp (double lnl_hom, double lnl_het) {
    return call_snp(2.0 * (lnl_hom - lnl_het));
}

inline
double lnl_multinomial_model_hom (double total, double n1) {
    if (n1 == total)
        return 0.0;
    else if (n1 < 0.25 * total)
        return total * log(0.25); // With epsilon estimate bounded at 1.0
    else
        return n1 * log(n1/total) + (total-n1) * log( (total-n1)/(3.0*total) );
}

inline
double lnl_multinomial_model_het (double total, double n1n2) {
    if (n1n2 == total)
        return total * log(0.5);
    else if (n1n2 < 0.5 * total)
        return total * log(0.25); // With epsilon estimate bounded at 1.0
    else
        return n1n2 * log( n1n2/(2.0*total) ) + (total-n1n2) * log( (total-n1n2)/(2.0*total) );
}

inline
double lr_multinomial_model_legacy (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {
    //
    // This function is to check that the refactored function gives the same
    // results as the original code (i.e. this code).
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);

    double l_ratio = (nuc_1 * log(nuc_1 / total));

    if (total - nuc_1 > 0.0)
        l_ratio += ((total - nuc_1) * log((total - nuc_1) / (3.0 * total)));

    if (nuc_1 + nuc_2 > 0.0)
        l_ratio -= ((nuc_1 + nuc_2) * log((nuc_1 + nuc_2) / (2.0 * total)));

    if (nuc_3 + nuc_4 > 0.0)
        l_ratio -= ((nuc_3 + nuc_4) * log((nuc_3 + nuc_4) / (2.0 * total)));

    l_ratio *= 2.0;

    return l_ratio;
}

inline
double lr_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {
    //
    // Method of Paul Hohenlohe <hohenlohe@uidaho.edu>, personal communication.
    //
    // For a diploid individual, there are ten possible genotypes
    // (four homozygous and six heterozygous genotypes).  We calculate
    // the likelihood of each possible genotype by using a multinomial
    // sampling distribution, which gives the probability of observing
    // a set of read counts (n1,n2,n3,n4) given a particular genotype.
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);
    assert(nuc_1 >= nuc_2 && nuc_2 >= nuc_3 && nuc_3 >= nuc_4);

    double l_ratio = 2.0 * (lnl_multinomial_model_hom(total, nuc_1) - lnl_multinomial_model_het(total, nuc_1+nuc_2));

    #ifdef DEBUG
    double l_ratio_legacy = lr_multinomial_model_legacy(nuc_1,nuc_2,nuc_3,nuc_4);
    assert( (l_ratio == 0.0 && l_ratio_legacy == 0.0) || almost_equal(l_ratio, l_ratio_legacy));
    #endif
    return l_ratio;
}

inline
double lr_bounded_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {

    //
    // Method of Paul Hohenlohe <hohenlohe@uidaho.edu>, personal communication.
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);

    //
    // Calculate the site specific error rate for homozygous and heterozygous genotypes.
    //
    double epsilon_hom  = (4.0 / 3.0) * ((total - nuc_1) / total);
    double epsilon_het  = 2.0 * ((nuc_3 + nuc_4) / total);

    //
    // Check if the error rate is above or below the specified bound.
    //
    if (epsilon_hom < bound_low)
        epsilon_hom = bound_low;
    else if (epsilon_hom > bound_high)
        epsilon_hom = bound_high;

    if (epsilon_het < bound_low)
        epsilon_het = bound_low;
    else if (epsilon_het > bound_high)
        epsilon_het = bound_high;

    //
    // Calculate the log likelihood for the homozygous and heterozygous genotypes.
    //
    double ln_L_hom = nuc_1 * log(1 - ((3.0/4.0) * epsilon_hom));
    ln_L_hom += epsilon_hom > 0.0 ? ((nuc_2 + nuc_3 + nuc_4) * log(epsilon_hom / 4.0)) : 0.0;

    double ln_L_het = (nuc_1 + nuc_2) * log(0.5 - (epsilon_het / 4.0));
    ln_L_het += epsilon_het > 0.0 ? ((nuc_3 + nuc_4) * log(epsilon_het / 4.0)) : 0.0;

    //
    // Calculate the likelihood ratio.
    //
    double l_ratio  = 2.0 * (ln_L_hom - ln_L_het);

    // cerr << "  Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Nuc_3: " << nuc_3 << " Nuc_4: " << nuc_4
    //   << " epsilon homozygote: " << epsilon_hom
    //   << " epsilon heterozygote: " << epsilon_het
    //   << " Log likelihood hom: " << ln_L_hom
    //   << " Log likelihood het: " << ln_L_het
    //   << " Likelihood ratio: " << l_ratio << "\n";

    return l_ratio;
}

inline
void SampleCall::set_call(snp_type c, Nt2 rank0_nt, Nt2 rank1_nt, long gq) {
    assert(gq >= 0 && gq <= 40);
    call_ = c;
    if (call_ == snp_type_hom) {
        nts_[0] = rank0_nt;
        gq_ = gq;
    } else if (call_ == snp_type_het) {
        if (rank0_nt < rank1_nt)
            nts_ = {{rank0_nt, rank1_nt}};
        else
            nts_ = {{rank1_nt, rank0_nt}};
        gq_ = gq;
    }
}

#endif // __MODELS_H__
