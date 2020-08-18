// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010 - 2012, Julian Catchen <jcatchen@uoregon.edu>
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

//
// models.cc -- routines to detect polymorphism (snp) and detect a lack of polymorphism (fixed).
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "models.h"

using namespace std;

vector<pair<char, int>> sort_acgt(const map<char, int>&);
void record_snp(SNP& snp, snp_type type, uint col, double l_ratio, const vector<pair<char, int>>& nuc);
void record_dummy_snp(SNP& snp, uint col);

int    barcode_size     = 5;
double barcode_err_freq = 0.0;

double heterozygote_limit = -qchisq(0.05,1);
double homozygote_limit   = qchisq(0.05,1);
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;

const map<string,modelt> model_strings = {
    {"snp", ::snp},
    {"bounded", ::bounded},
    {"fixed", ::fixed},
    {"marukihigh", ::marukihigh},
    {"marukilow", ::marukilow},
};

// Quantiles of the Chi2 distribution for VCF's GQ.
//$ Rscript -e "cat(paste(qchisq(1-10^(-(1:40)/10), 1), '\n', sep=''))"
const std::array<double,40> chisq1ddf_gq = {{
    0.0679615365255401, 0.230764766661193, 0.452421556097698, 0.714036104152315, 1.00448471501902,
    1.31668063342863, 1.64583886397469, 1.98858107774449, 2.34243609502603, 2.70554345409542,
    3.07646857913928, 3.45408274293754, 3.83748240184564, 4.22593339019369, 4.61883133337202,
    5.01567294625387, 5.41603482049954, 5.81955747770787, 6.22593319775977, 6.63489660102121,
    7.04621727098416, 7.45969391025114, 7.87514966368955, 8.29242834051146, 8.71139133617301,
    9.13191510451069, 9.55388906647803, 9.97721386826398, 10.4017999212068, 10.8275661706627,
    11.2544390521783, 11.6823516018745, 12.1112426945654, 12.5410563882771, 12.9717413578738,
    13.403250403671, 13.8355400234795, 14.2685700385079, 14.7023032652319, 15.1367052266236
}};

double qchisq(double alpha, size_t df) {
    //
    // Quantiles (1-alpha) of the Chi2 distribution, as given by
    //$ alphas="1e-1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6"
    //$ df=1
    //$ R -e "for(a in c($alphas)) {q=qchisq(1-a, $df); cat(paste('{',a,',',q,'}, ',sep=''));}"
    //
    // The value for a=0.05,df=1 is 3.84 for backward compatiblity (instead of
    // 3.84145882069412).
    //
    static const vector<map<double,double>> table = {
        { //df=1
            {0.1,2.70554345409542}, {0.05,3.84}, {0.01,6.63489660102121},
            {0.005,7.87943857662241}, {0.001,10.8275661706627}, {0.0005,12.1156651463974},
            {0.0001,15.1367052266236}, {0.00005,16.4481102100082}, {0.00001,19.5114209646663},
            {0.000005,20.8372870225104}, {0.000001,23.9281269768795}
        },{ //df=2
            {0.1,4.60517018598809}, {0.05,5.99146454710798}, {0.01,9.21034037197618},
            {0.005,10.5966347330961}, {0.001,13.8155105579643}, {0.0005,15.2018049190844},
            {0.0001,18.4206807439526}, {0.00005,19.8069751050725}, {0.00001,23.0258509299496},
            {0.000005,24.4121452910472}, {0.000001,27.631021115871}
        }
    };

    assert(df > 0 && df <= table.size());
    try {
        return table[df-1].at(alpha);
    } catch (std::out_of_range&) {
        cerr << "Error: Unsupported alpha value '" << alpha << "'; pick one of 0.05, 0.01, 0.005, 0.001...\n";
        throw std::invalid_argument("qchisq");
    }
}

modelt parse_model_type(const string& arg) {
    try {
        return model_strings.at(arg);
    } catch (std::out_of_range&) {
        cerr << "Error: Unknown model '" << arg << "'.\n";
        throw std::invalid_argument("set_model_type");
    }
}

void set_model_thresholds(double alpha) {
    double x = qchisq(alpha, 1);
    homozygote_limit = x;
    heterozygote_limit = -x;
}

string to_string(modelt model_type) {
    for (auto& m : model_strings)
        if (model_type == m.second)
            return m.first;
    DOES_NOT_HAPPEN;
    return string();
}

vector<pair<char, int>> sort_acgt(const map<char, int>& counts) {
    vector<pair<char, int>> nuc;
    for (map<char, int>::const_iterator i=counts.begin(); i!=counts.end(); i++)
        if (i->first != 'N')
            nuc.push_back(make_pair(i->first, i->second));
    sort(nuc.begin(), nuc.end(), compare_pair);
    return nuc;
}

void record_snp(SNP& snp, snp_type type, uint col, double l_ratio, const vector<pair<char, int>>& nuc) {

    snp.type   = type;
    snp.col    = col;
    snp.lratio = l_ratio;
    snp.rank_1 = nuc[0].first;

    switch (type) {
    case snp_type_het:
        snp.rank_2 = nuc[1].first;
        break;
    case snp_type_hom:
        snp.rank_2 = '-';
        break;
    default:
        // snp_type_unk, with at least one observation (otherwise record_dummy_snp
        // would have been called)
        // A rank_2 nucleotide is set only if at least two nucleotides were observed.
        snp.rank_2 = nuc[1].second > 0 ? nuc[1].first : '-';
        break;
    }

    snp.rank_3 = 0;
    snp.rank_4 = 0;
}

void record_dummy_snp(SNP& snp, uint col) {
    snp.type   = snp_type_unk;
    snp.col    = col;
    snp.lratio = 0.0;
    snp.rank_1 = 'N';
    snp.rank_2 = '-';
    snp.rank_3 = 0;
    snp.rank_4 = 0;
}

void call_multinomial_snp(MergedStack *tag, int col, map<char, int> &n, bool record_snps) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        if (record_snps) {
            tag->snps.push_back(new SNP());
            record_dummy_snp(*tag->snps.back(), col);
        }
        return;
    }

    double l_ratio = lr_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    if (record_snps) {
        tag->snps.push_back(new SNP());
        record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
    }
}

void call_bounded_multinomial_snp(MergedStack *tag, int col, map<char, int> &n, bool record_snps) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        if (record_snps) {
            tag->snps.push_back(new SNP());
            record_dummy_snp(*tag->snps.back(), col);
        }
        return;
    }

    double l_ratio = lr_bounded_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    if (record_snps) {
        tag->snps.push_back(new SNP());
        record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
    }
}

void call_multinomial_snp(Locus *tag, int col, map<char, int> &n) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        record_dummy_snp(*tag->snps[col], col);
        return;
    }

    double l_ratio = lr_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    record_snp(*tag->snps[col], type, col, l_ratio, nuc);
}

void call_bounded_multinomial_snp(Locus *tag, int col, map<char, int> &n) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        record_dummy_snp(*tag->snps[col], col);
        return;
    }

    double l_ratio = lr_bounded_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    record_snp(*tag->snps[col], type, col, l_ratio, nuc);
}

void call_multinomial_fixed (MergedStack *tag, int col, map<char, int> &n) {
    const double nucleotide_fixed_limit = 1.92;

    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        tag->snps.push_back(new SNP());
        record_dummy_snp(*tag->snps.back(), col);
        return;
    }

    double l_ratio;
    double nuc_1 = nuc[0].second;
    double nuc_2 = nuc[1].second;
    {
        //
        // Method of Paul Hohenlohe <hohenlo@uoregon.edu>, personal communication.
        //
        // Each population sample contains DNA from 6 individuals, so a
        // sample of 12 alleles from the population. We want to assign a
        // nucleotide (A,C,G,T) to each position where the population is
        // fixed or nearly so, and N to each position that is either
        // polymorphic within the population or has insufficient coverage
        // depth to make a call. We can do this with a likelihood ratio
        // test of the read counts, testing whether the allele frequency
        // of the dominant allele is significantly larger than some
        // threshold p) , stepping through each nucleotide position across
        // RAD tags.
        //

        double epsilon = -1 * (log(1 - barcode_err_freq) / barcode_size);

        l_ratio  =
            nuc_1 * log( ((4 * nuc_1 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) /
                         ((4 * p_freq * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

        l_ratio +=
            nuc_2 * log( ((4 * nuc_2 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) /
                         ((4 * (1 - p_freq) * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

        //cerr << "Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Likelihood ratio: " << l_ratio << "\n";
    }

    double n_ratio = nuc_1 / (nuc_1 + nuc_2);

    snp_type type = n_ratio < p_freq || l_ratio < nucleotide_fixed_limit ? snp_type_unk : snp_type_hom;

    tag->snps.push_back(new SNP());
    record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
}

//
// ln L(1/2) = ln(n! / n_1!n_2!n_3!n_4!) +
//               (n_1 + n_2) * ln(n_1 + n_2 / 2n) +
//               (n_3 + n_4) * ln(n_3 + n_4 / 2n)
//
double
heterozygous_likelihood(int col, map<char, int> &nuc)
{
    vector<pair<char, int> > cnts;
    map<char, int>::iterator i;

    double n = 0;
    for (i = nuc.begin(); i != nuc.end(); i++) {
        n += i->second;
        cnts.push_back(make_pair(i->first, i->second));
    }

    sort(cnts.begin(), cnts.end(), compare_pair);

    double n_1 = cnts[0].second;
    double n_2 = cnts[1].second;
    double n_3 = cnts[2].second;
    double n_4 = cnts[3].second;

    double term_1 =
        reduced_log_factorial(n, n_1) -
        (log_factorial(n_2) + log_factorial(n_3) + log_factorial(n_4));

    double term_3 = (n_3 + n_4 > 0) ? log((n_3 + n_4) / (2 * n)) : 0;

    double lnl =
        term_1 +
        ((n_1 + n_2) * log((n_1 + n_2) / (2 * n))) +
        ((n_3 + n_4) * term_3);

    return lnl;
}

//
// ln L(1/1) = ln(n! / n_1!n_2!n_3!n_4!) +
//               n_1 * ln(n_1 / n) +
//               (n - n_1) * ln(n - n_1 / 3n)
//
double
homozygous_likelihood(int col, map<char, int> &nuc)
{
    vector<pair<char, int> > cnts;
    map<char, int>::iterator i;

    double n = 0;
    for (i = nuc.begin(); i != nuc.end(); i++) {
        n += i->second;
        cnts.push_back(make_pair(i->first, i->second));
    }

    sort(cnts.begin(), cnts.end(), compare_pair);

    double n_1 = cnts[0].second;
    double n_2 = cnts[1].second;
    double n_3 = cnts[2].second;
    double n_4 = cnts[3].second;

    double term_1 =
        reduced_log_factorial(n, n_1) -
        (log_factorial(n_2) + log_factorial(n_3) + log_factorial(n_4));

    double term_3 = n - n_1 > 0 ? log((n - n_1) / (3 * n)) : 0;

    double lnl =
        term_1 +
        (n_1 * log(n_1 / n)) +
        ((n - n_1) * term_3);

    return lnl;
}

Nt2 SiteCall::most_frequent_allele() const {
    assert(!alleles().empty());
    auto a = alleles().begin();
    auto best = a;
    for(; a!=alleles().end(); ++a)
        if (a->second > best->second)
            best = a;
    return best->first;
}

void SiteCall::filter_mac(size_t min_mac) {
    assert(alleles().size() > 1);
    Counts<Nt2> allele_counts = tally_allele_counts(sample_calls());
    if (allele_counts.sorted()[1].first < min_mac)
        *this = SiteCall(most_frequent_allele());
}

Counts<Nt2> SiteCall::tally_allele_counts(const vector<SampleCall>& spldata) {
    Counts<Nt2> counts;
    for (const SampleCall& sd : spldata) {
        switch (sd.call()) {
        case snp_type_hom :
            counts.increment(sd.nt0());
            counts.increment(sd.nt0());
            break;
        case snp_type_het :
            counts.increment(sd.nt0());
            counts.increment(sd.nt1());
            break;
        default:
            // snp_type_unk
            break;
        }
    }
    return counts;
}

map<Nt2,double> SiteCall::tally_allele_freqs(const vector<SampleCall>& spldata) {
    //
    // Tally the existing alleles and their frequencies.
    //
    map<Nt2,double> allele_freqs;
    Counts<Nt2> counts = tally_allele_counts(spldata);
    array<pair<size_t,Nt2>,4> sorted_alleles = counts.sorted();
    size_t tot = counts.sum();
    for (auto& a : sorted_alleles) {
        if (a.first == 0)
            break;
        allele_freqs.insert({a.second, double(a.first)/tot});
    }
    return allele_freqs;
}

ostream& operator<<(ostream& os, const SampleCall& c) {
    ostream copy (os.rdbuf());
    copy << std::setprecision(4);
    // Genotype.
    switch(c.call()) {
    case snp_type_hom: os << c.nt0() << "/" << c.nt0(); break;
    case snp_type_het: os << std::min(c.nt0(), c.nt1()) << "/" << std::max(c.nt0(), c.nt1()); break;
    case snp_type_unk: os << "u"; break;
    case snp_type_discarded: os << "d"; break;
    default: DOES_NOT_HAPPEN; break;
    }
    // Likelihoods.
    os << "\t{" << c.lnls() << "}";
    return os;
}

void SiteCall::print(ostream& os, const SiteCounts& depths) {
    // Total depths.
    os << "tot depths {" << depths.tot << "}\n";
    // Alleles.
    os << "alleles";
    for (auto& a : alleles())
        os << " " << a.first << "(" << a.second << ")";
    os << "\n";
    // Samples.
    os << "samples";
    for (size_t s=0; s<depths.samples.size(); ++s) {
        // Index.
        os << "\n" << s;
        auto& dp = depths.samples[s];
        if (dp.sum() == 0) {
            os << "\t.";
        } else {
            // Depths.
            os << "\t{" << dp << "}";
            if (alleles().size() >= 2)
                os << "\t" << sample_calls()[s];
        }
    }
}

SiteCall MultinomialModel::call(const SiteCounts& depths) const {

    size_t n_samples = depths.mpopi->samples().size();

    //
    // Make genotype calls.
    //
    vector<SampleCall> sample_calls (n_samples);
    array<pair<size_t,Nt2>,4> sorted;
    for (size_t sample=0; sample<n_samples; ++sample) {
        const Counts<Nt2>& sdepths = depths.samples[sample];
        size_t dp = sdepths.sum();
        if (dp == 0)
            continue;

        sorted = sdepths.sorted();
        SampleCall& c = sample_calls[sample];
        double lnl_hom = lnl_multinomial_model_hom(dp, sorted[0].first);
        double lnl_het = lnl_multinomial_model_het(dp, sorted[0].first+sorted[1].first);
        c.lnls().set(sorted[0].second, sorted[0].second, lnl_hom);
        c.lnls().set(sorted[0].second, sorted[1].second, lnl_het);

        snp_type gt_call = call_snp(lnl_hom, lnl_het);
        long gq = vcf_lnl2gq(std::max(lnl_hom, lnl_het), std::min(lnl_hom, lnl_het));
        c.set_call(gt_call, sorted[0].second, sorted[1].second, gq);
    }

    //
    // Record the existing alleles and their frequencies.
    //
    map<Nt2,double> allele_freqs = SiteCall::tally_allele_freqs(sample_calls);

    //
    // Compute the missing genotype likelihoods, if any.
    //
    if (allele_freqs.size() < 2) {
        sample_calls = vector<SampleCall>();
    } else {
        for (size_t sample=0; sample<n_samples; ++sample) {
            const Counts<Nt2>& sdepths = depths.samples[sample];
            size_t dp = sdepths.sum();
            if (dp == 0)
                continue;

            SampleCall& c = sample_calls[sample];
            for (auto nt1_it=allele_freqs.begin(); nt1_it!=allele_freqs.end(); ++nt1_it) {
                Nt2 nt1 (nt1_it->first);
                // Homozygote.
                if (!c.lnls().has_lik(nt1, nt1)) {
                    double lnl = lnl_multinomial_model_hom(dp, sdepths[nt1]);
                    c.lnls().set(nt1, nt1, lnl);
                }
                // Heterozygote(s).
                auto nt2_it = nt1_it;
                ++nt2_it;
                for (; nt2_it!=allele_freqs.end(); ++nt2_it) {
                    Nt2 nt2 (nt2_it->first);
                    if (!c.lnls().has_lik(nt1, nt2)) {
                        double lnl = lnl_multinomial_model_het(dp, sdepths[nt1]+sdepths[nt2]);
                        c.lnls().set(nt1, nt2, lnl);
                    }
                }
            }
        }
    }

    return SiteCall(move(allele_freqs), move(sample_calls));
}

double MarukiHighModel::calc_hom_lnl(double n, double n1) const {
    // This returns the same value as the Hohenlohe ('snp/binomial') model except
    // when the error rate estimate is bounded at 1.0 (i.e. `n1 < 0.25*n` c.f.
    // Hohenlohe equations).
    if (n1 == n)
        return 0.0;
    else if (n1 == 0.0)
        return n * log(1.0/3.0);
    else
        return n1 * log(n1/n) + (n-n1) * log( (n-n1)/(3.0*n) );
}

double MarukiHighModel::calc_het_lnl(double n, double n1n2) const {
    // This returns the same value as the Hohenlohe ('snp/binomial') model except
    // when the error rate estimate is bounded at 1.0 (i.e. `n1n2 < 0.5*n` c.f.
    // Hohenlohe equations).
    if (n1n2 == n)
        return n * log(0.5);
    else if (n1n2 < (1.0/3.0) * n)
        return n1n2 * log(1.0/6.0) + (n-n1n2) * log(1.0/3.0);
    else
        return n1n2 * log( n1n2/(2.0*n) ) + (n-n1n2) * log((n-n1n2)/(2.0*n) );
}

SiteCall MarukiHighModel::call(const SiteCounts& depths) const {

    /*
     * For this model the procedure is:
     * I. Obtain the most commonly seen nucleotide M.
     * II. Look for alternative alleles, if any:
     *     For each sample:
     *         If there is a genotype significantly better than MM:
     *             (This genotype can be Mm, mm or mn.)
     *             The site is polymorphic.
     *             Record m as an alternative allele.
     *             If the genotype is mn AND is significantly better than mm:
     *                 Record n as an allele.
     * III. Given the known alleles, compute the likelihoods for all possible
     *      genotypes (n.b. most of the non-trivial ones have already been
     *      computed), and call genotypes.
     */

    const size_t n_samples = depths.mpopi->samples().size();

    if (depths.tot.sum() == 0)
        return SiteCall(map<Nt2,double>(), vector<SampleCall>());

    //
    // I.
    // Count the observed nucleotides of the site for all samples; set
    // `SampleCall::depths_`.
    // Then find the most common nucleotide across the population.
    //
    Nt2 nt_ref = depths.tot.sorted()[0].second;

    //
    // II.
    // Look for alternative alleles; start filling SampleCall::lnls_.
    //
    set<Nt2> alleles;
    vector<SampleCall> sample_calls (n_samples);
    alleles.insert(nt_ref);
    array<pair<size_t,Nt2>,4> sorted;
    for (size_t sample=0; sample<n_samples; ++sample) {
        const Counts<Nt2>& sdepths = depths.samples[sample];
        size_t dp = sdepths.sum();
        if (dp == 0)
            continue;

        // Find the best genotype for the sample -- this is either the homzygote
        // for the rank0 nucleotide or the heterozygote for the rank0 and rank1
        // nucleotides.
        sorted = sdepths.sorted();
        Nt2 nt0 = sorted[0].second;
        Nt2 nt1 = sorted[1].second;
        size_t dp0 = sorted[0].first;
        size_t dp1 = sorted[1].first;

        if (dp1 == 0 && nt0 == nt_ref)
            // We can wait until we know whether the site is fixed.
            continue;

        SampleCall& c = sample_calls[sample];
        double lnl_hom = calc_hom_lnl(dp, dp0);
        double lnl_het = calc_het_lnl(dp, dp0+dp1);
        c.lnls().set(nt0, nt0, lnl_hom);
        c.lnls().set(nt0, nt1, lnl_het);

        // Make sure the sample would have a significant genotype call provided
        // the site was polymorphic (otherwise a genotype can be significant in
        // comparison with the ref homozygote but not in comparison with the
        // second best genotype). This is a slight modification to Maruki &
        // Lynch's method to avoid low-coverage weirdnesses. Note that the
        // polymorphism discovery alpha could be set lower than the genotype call
        // alpha (e.g. to account for multiple testing).
        if (lnl_hom > lnl_het) {
            if(!lrtest(lnl_hom, lnl_het, gt_threshold_))
                continue;
        } else {
            if(!lrtest(lnl_het, lnl_hom, gt_threshold_))
                continue;
        }

        // Compare this likelihood to that of the ref,ref homozygote.
        if (nt0 == nt_ref) {
            if (lrtest(lnl_het, lnl_hom, var_threshold_))
                // Record the alternative allele.
                alleles.insert(nt1);
        } else {
            double lnl_ref = calc_hom_lnl(dp, sdepths[nt_ref]);
            c.lnls().set(nt_ref, nt_ref, lnl_ref);
            double lnl_best = std::max(lnl_hom, lnl_het);
            if (lrtest(lnl_best, lnl_ref, var_threshold_)) {
                // Record one alternative allele.
                alleles.insert(nt0);
                if (nt1!=nt_ref && lrtest(lnl_het, lnl_hom, var_threshold_))
                    // Record a second alternative allele (the SNP is at least ternary).
                    alleles.insert(nt1);
            }
        }
    }

    //
    // III.
    // Compute the likelihoods for all possible genotypes & call genotypes.
    //
    if (alleles.size() == 1) {
        sample_calls = vector<SampleCall>();
    } else {
        for (size_t sample=0; sample<n_samples; ++sample) {
            const Counts<Nt2>& sdepths = depths.samples[sample];
            size_t dp = sdepths.sum();
            if (dp == 0)
                continue;

            SampleCall& c = sample_calls[sample];
            for (auto nt1=alleles.begin(); nt1!=alleles.end(); ++nt1) {
                // Homozygote.
                if (!c.lnls().has_lik(*nt1, *nt1))
                    c.lnls().set(*nt1, *nt1, calc_hom_lnl(dp, sdepths[*nt1]));
                // Heterozygote(s).
                auto nt2 = nt1;
                ++nt2;
                for (; nt2!=alleles.end(); ++nt2)
                    if (!c.lnls().has_lik(*nt1, *nt2))
                        c.lnls().set(*nt1, *nt2, calc_het_lnl(dp, sdepths[*nt1]+sdepths[*nt2]));
            }

            // Call the genotype -- skiping ignored alleles.
            sorted = sdepths.sorted();
            auto nt0 = sorted.begin();
            while(!alleles.count(nt0->second)) {
                ++nt0;
                assert(nt0 != sorted.end());
            }
            auto nt1 = nt0;
            ++nt1;
            while(!alleles.count(nt1->second)) {
                ++nt1;
                assert(nt1 != sorted.end());
            }
            double lnl_hom = c.lnls().at(nt0->second, nt0->second);
            double lnl_het = c.lnls().at(nt0->second, nt1->second);
            snp_type call;
            long gq;
            if (lnl_hom > lnl_het) {
                call = lrtest(lnl_hom, lnl_het, gt_threshold_) ? snp_type_hom : snp_type_unk;
                gq = vcf_lnl2gq(lnl_hom, lnl_het);
            } else {
                call = lrtest(lnl_het, lnl_hom, gt_threshold_) ? snp_type_het : snp_type_unk;
                gq = vcf_lnl2gq(lnl_het, lnl_hom);
            }

            c.set_call(call, nt0->second, nt1->second, gq);
        }
    }

    //
    // Finally, compute allele frequencies.
    //
    map<Nt2,double> allele_freqs;
    if (alleles.size() == 1) {
        allele_freqs.insert({*alleles.begin(), 1.0});
    } else {
        allele_freqs = SiteCall::tally_allele_freqs(sample_calls);
        assert(allele_freqs.size() == alleles.size()
               || (allele_freqs.size() == alleles.size()-1 && !allele_freqs.count(nt_ref)));
        // Note: There are limit cases (esp. low coverage) where nt_ref does
        // not appear in any of the significant genotypes.
    }

    return SiteCall(move(allele_freqs), move(sample_calls));
}

double MarukiLowModel::calc_fixed_lnl(double n_tot, double n_M_tot) const {
    assert(n_tot > 0.0);
    assert(n_M_tot > 0.0);
    if (n_M_tot == n_tot)
        return 0.0;
    else
        return n_M_tot * log(n_M_tot/n_tot) + (n_tot-n_M_tot) * log((n_tot-n_M_tot)/n_tot/3.0);
}

double MarukiLowModel::calc_dimorph_lnl(double freq_MM, double freq_Mm, double freq_mm, const vector<LikData>& liks) const {
    double lnl = 0.0;
    // Sum over samples.
    for (const LikData& s_liks : liks)
        if (s_liks.has_data)
            // If !has_data, the sum is 1 and its log 0.
            lnl += calc_ln_weighted_sum(freq_MM, freq_Mm, freq_mm, s_liks);
    assert(lnl <= 0.0);
    return lnl;
}

double MarukiLowModel::calc_ln_weighted_sum(double freq_MM, double freq_Mm, double freq_mm, const LikData& s_liks) const {
    double weighted_sum = freq_MM * s_liks.l_MM + freq_Mm * s_liks.l_Mm + freq_mm * s_liks.l_mm;
    #pragma omp atomic
    ++n_wsum_tot_;
    if (weighted_sum >= std::numeric_limits<double>::min()) {
        return log(weighted_sum);
    } else {
        // `weigted_sum` is subnormal or zero.
        #pragma omp atomic
        ++n_wsum_underflows_;
        return calc_ln_weighted_sum_safe(freq_MM, freq_Mm, freq_mm, s_liks);
    }
}

double MarukiLowModel::calc_ln_weighted_sum_safe(double freq_MM, double freq_Mm, double freq_mm, const LikData& s_liks) const {
    array<pair<double,double>,3> s = {{
        {s_liks.lnl_MM, freq_MM},
        {s_liks.lnl_Mm, freq_Mm},
        {s_liks.lnl_mm, freq_mm}
    }};
    std::sort(s.begin(), s.end()); // `s` is sorted by increasing lnl.
    if (s[2].second > 0.0)
        return s[2].first + log(s[2].second
                                + s[1].second * exp(s[1].first-s[2].first)
                                + s[0].second * exp(s[0].first-s[2].first)
                                );
    else if (s[1].second > 0.0)
        return s[1].first + log(s[1].second + s[0].second * exp(s[0].first-s[1].first));
    else
        return s[0].first + log(s[0].second);
}

SiteCall MarukiLowModel::call(const SiteCounts& depths) const {

    /*
     * For this model the procedure is:
     * I. Compute the maximum likelihood for the fixed-site hypothesis (straightforward).
     * II. Compute the maximum likelihood for the dimorphic-site hypothesis and
     *     estimate the genotype frequencies:
     *     1. Compute the error rate estimate.
     *     2. Compute the base likelihoods for all samples, all genotypes.
     *     3. Compute the starting major allele frequency (`p`) value.
     *     4. Compute the likelihoods for the site for all possible disequilibrium
     *         coefficients (`d_a`); keep the best one.
     *     5. Optimize the genotype frequencies by finding a local maximum.
     * III. Test whether the site is polymorphic.
     * IV. Compute the likelihoods for all samples, all genotypes using Bayes'
     *      theorem, and call genotypes.
     */

    const size_t n_samples = depths.mpopi->samples().size();

    size_t dp_tot = depths.tot.sum();
    if (dp_tot == 0)
        return SiteCall();

    //
    // I. Likelihood for the fixed-site hypothesis.
    //

    array<pair<size_t,Nt2>,4> sorted = depths.tot.sorted();
    Nt2 nt_M = sorted[0].second;
    Nt2 nt_m = sorted[1].second;
    size_t n_M_tot = sorted[0].first;
    size_t n_m_tot = sorted[1].first;

    double lnl_fixed = calc_fixed_lnl(dp_tot, n_M_tot);

    if (!lrtest(0.0, lnl_fixed, var_threshold_))
        // Fixed: `lnl_fixed` is high enough than even when `lnl_dimorphic` is
        // at its maximum value of 0, dimorphism isn't significant. This happens
        // when there are either very few non-major-allele reads or very few reads
        // overall.
        return SiteCall(nt_M);

    //
    // II. Compute and optimize the likelihood for the dimorphic-site hypothesis
    //

    // 1. Error rate.
    // Note: @Nick, Apr 2017: Maruki's description uses a simple estimator where
    // the number of invisible errors (from M to m and inversely) is estimated
    // as being half the number of observable errors. But this when few
    // observable errors this is an underestimate as ideally we would like to
    // integrate; in particular the error rate is estimated to zero when
    // there aren't any observable errors and this leads to meaningless null
    // likelihoods. Doing the integration isn't practical but we can reduce
    // the problem by adding (see edit) to the numerator and denominator.
    // Edit. Oct 2017: 1 is way too large; adding 0.1 yields a more reasonable
    // prior.
    double e = 3.0 / 2.0 * (dp_tot - n_M_tot - n_m_tot + 0.1) / double(dp_tot + 0.1);
    assert(e > 0.0 && e < 1.0);
    #pragma omp atomic
    ++n_called_sites_;
    #pragma omp atomic
    sum_site_err_rates_ += e;

    // 2. Base likelihoods.
    vector<LikData> liks;
    {
        liks.reserve(n_samples);
        double ln_err_hom = log(e/3.0);
        double ln_hit_hom = log(1-e);
        double ln_err_het = log(e/3.0);
        double ln_hit_het = log(0.5-e/3.0);
        for (size_t sample=0; sample<n_samples; ++sample) {
            const Counts<Nt2>& sdepths = depths.samples[sample];
            size_t dp = sdepths.sum();
            if (dp == 0) {
                liks.push_back(LikData());
                continue;
            }

            double lnl_MM = sdepths[nt_M] * ln_hit_hom + (dp-sdepths[nt_M]) * ln_err_hom;
            double lnl_mm = sdepths[nt_m] * ln_hit_hom + (dp-sdepths[nt_m]) * ln_err_hom;
            double lnl_Mm = (sdepths[nt_M]+sdepths[nt_m]) * ln_hit_het + (dp-(sdepths[nt_M]+sdepths[nt_m])) * ln_err_het;;
            liks.push_back(LikData(lnl_MM, lnl_Mm, lnl_mm));
        }
    }

    // 3. Initial major allele frequency.
    double p = (2*n_M_tot + 2*n_m_tot > dp_tot) ?
            0.5 * double(3*n_M_tot + n_m_tot - dp_tot) / (2*n_M_tot + 2*n_m_tot - dp_tot)
            : 0.5; // if `n_M_tot == n_m_tot == dp_tot / 4`.

    // 4. Initial disequilibrium coefficient.
    auto flush_freq = [](double& f) {
        if (f < 1e-12)
            f = 0.0;
        else if (f > 1.0-1e-12)
            f = 1.0;
    };
    double d_a_best = 0.0;
    double lnl_dimorph = std::numeric_limits<double>::lowest();
    {
        double d_a_min = max(-p*p, -(1-p)*(1-p));
        double d_a_max = p*(1-p);
        for (double d_a = d_a_min; d_a < d_a_max + 1e-12; d_a += 1.0/n_samples) {
            double f_MM = p*p + d_a;
            double f_Mm = 2 * (p*(1-p) - d_a);
            double f_mm = (1-p) * (1-p) + d_a;
            flush_freq(f_MM);
            flush_freq(f_Mm);
            flush_freq(f_mm);
            double lnl = calc_dimorph_lnl(f_MM, f_Mm, f_mm, liks);
            if (lnl > lnl_dimorph) {
                d_a_best = d_a;
                lnl_dimorph = lnl;
            }
        }
    }
    assert(lnl_dimorph > std::numeric_limits<double>::lowest());

    // 5. Find a local maximum.
    double freq_MM = p*p + d_a_best;
    double freq_Mm = 2 * (p*(1-p) - d_a_best);
    double freq_mm = (1-p) * (1-p) + d_a_best;
    flush_freq(freq_MM);
    flush_freq(freq_Mm);
    flush_freq(freq_mm);
    double x = 1.0/n_samples;
    {
        double lnl_prev;
        do {
            lnl_prev = lnl_dimorph;
            array<array<double,3>,6> neighbors = {{
                {{freq_MM+x, freq_Mm-x, freq_mm}},
                {{freq_MM-x, freq_Mm+x, freq_mm}},
                {{freq_MM+x, freq_Mm,   freq_mm-x}},
                {{freq_MM-x, freq_Mm,   freq_mm+x}},
                {{freq_MM,   freq_Mm+x, freq_mm-x}},
                {{freq_MM,   freq_Mm-x, freq_mm+x}}
            }};
            for(auto& n : neighbors) {
                double f_MM = n[0];
                double f_Mm = n[1];
                double f_mm = n[2];
                if (f_MM < 0.0-1e-12 || f_MM > 1.0+1e-12
                 || f_Mm < 0.0-1e-12 || f_Mm > 1.0+1e-12
                 || f_mm < 0.0-1e-12 || f_mm > 1.0+1e-12)
                    // Out of bounds.
                    continue;
                flush_freq(f_MM);
                flush_freq(f_Mm);
                flush_freq(f_mm);

                double lnl = calc_dimorph_lnl(f_MM, f_Mm, f_mm, liks);
                if (lnl > lnl_dimorph) {
                    // Update the frequencies.
                    freq_MM = f_MM;
                    freq_Mm = f_Mm;
                    freq_mm = f_mm;
                    lnl_dimorph = lnl;
                }
            }
        } while (lnl_dimorph > lnl_prev);
    }
    assert(almost_equal(freq_MM+freq_Mm+freq_mm, 1.0));

    //
    // III. Test whether the site is polymorphic.
    //

    if (!lrtest(lnl_dimorph, lnl_fixed, var_threshold_))
        return SiteCall(nt_M);

    map<Nt2,double> allele_freqs = {
        {nt_M, freq_MM+0.5*freq_Mm},
        {nt_m, freq_mm+0.5*freq_Mm},
    };
    for (auto& a : allele_freqs)
        flush_freq(a.second);

    long dimorph_qual = vcf_lnl2gq(lnl_dimorph, lnl_fixed);

    //
    // IV. Bayes-corrected likelihoods & genotypes.
    //
    // Note: @Nick Oct 2017: Optimal genotype frequencies may be 0.0 when the MAC
    // or the number of samples is low. To avoid feeding null frequencies into
    // the genotype equations, which may lead to distorted likelihood ratios
    // (e.g. for a candidate het with more coverage for the minor allele, when
    // the minor homozygote frequency is 0.0), we add 1/n_samples to all three
    // frequencies.
    //

    vector<SampleCall> sample_calls (n_samples);
    {
        sample_calls.reserve(n_samples);
        freq_MM = (freq_MM + x) / (1.0 + 3 * x);
        freq_Mm = (freq_Mm + x) / (1.0 + 3 * x);
        freq_mm = (freq_mm + x) / (1.0 + 3 * x);
        double log_f_MM = log(freq_MM);
        double log_f_Mm = log(freq_Mm);
        double log_f_mm = log(freq_mm);
        size_t gt_MM = GtLiks::gt_index(nt_M, nt_M);
        size_t gt_Mm = GtLiks::gt_index(nt_M, nt_m);
        size_t gt_mm = GtLiks::gt_index(nt_m, nt_m);
        for (size_t sample=0; sample<n_samples; ++sample) {
            const LikData& s_liks = liks[sample];
            SampleCall& s_call = sample_calls[sample];

            double w_sum = calc_ln_weighted_sum(freq_MM, freq_Mm, freq_mm, s_liks);

            double lnl_MM = s_liks.lnl_MM + log_f_MM - w_sum;
            double lnl_Mm = s_liks.lnl_Mm + log_f_Mm - w_sum;
            double lnl_mm = s_liks.lnl_mm + log_f_mm - w_sum;
            assert(lnl_MM < 0.0+1e-12);
            assert(lnl_Mm < 0.0+1e-12);
            assert(lnl_mm < 0.0+1e-12);
            if (lnl_MM > 0.0)
                lnl_MM = 0.0;
            if (lnl_Mm > 0.0)
                lnl_Mm = 0.0;
            if (lnl_mm > 0.0)
                lnl_mm = 0.0;
            s_call.lnls().set(gt_MM, lnl_MM);
            s_call.lnls().set(gt_Mm, lnl_Mm);
            s_call.lnls().set(gt_mm, lnl_mm);

            array<pair<double,pair<Nt2,Nt2>>,3> lnls {{
                {lnl_MM, {nt_M,nt_M}},
                {lnl_Mm, {nt_M,nt_m}},
                {lnl_mm, {nt_m,nt_m}}
            }};
            sort(lnls.rbegin(), lnls.rend());

            if (lrtest(lnls[0].first, lnls[1].first, gt_threshold_)) {
                pair<Nt2,Nt2>& nts = lnls[0].second;
                long gq = vcf_lnl2gq(lnls[0].first, lnls[1].first);
                s_call.set_call(nts.first == nts.second ? snp_type_hom : snp_type_het, nts.first, nts.second, gq);
            }
        }
    }

    return SiteCall(move(allele_freqs), dimorph_qual, move(sample_calls));
}
