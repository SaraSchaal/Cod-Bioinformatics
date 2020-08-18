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
#include <algorithm>
#include <vector>
#include <iomanip>

#include "export_formats.h"

#include "utils.h"
#include "sql_utilities.h"
#include "MetaPopInfo.h"
#include "nucleotides.h"
#include "Vcf.h"
#include "ordered.h"

using namespace std;

extern string      out_path;
extern string      out_prefix;
extern bool        loci_ordered;
extern bool        ordered_export;

void tally_complete_haplotypes(
        Datum const*const* data,
        size_t n_samples,
        strand_type loc_strand, // (Alphabetic order is used for frequency ties.)
        vector<pair<const char*, size_t>>& haps_sorted_decr_freq,
        map<const char*, size_t, LessCStrs>& hap_indexes_map
) {
    haps_sorted_decr_freq.clear();
    hap_indexes_map.clear();

    // Tally the haplotype frequencies.
    // N.B. We first use the map to store freqs/counts; later we replace
    // those values with the indexes.
    for (size_t sample=0; sample<n_samples; ++sample) {
        if (data[sample] == NULL)
            continue;
        const vector<char*>& s_haps = data[sample]->obshap;
        assert(s_haps.size() == 2);
        if (strchr(s_haps[0], 'N') != NULL) {
            // Partially unresolved.
            assert(strchr(s_haps[1], 'N') != NULL);
            continue;
        }
        assert(strchr(s_haps[1], 'N') == NULL);
        ++hap_indexes_map[s_haps[0]];
        ++hap_indexes_map[s_haps[1]];
    }

    // Sort the haplotypes.
    haps_sorted_decr_freq.assign(hap_indexes_map.begin(), hap_indexes_map.end());
    std::sort(haps_sorted_decr_freq.begin(), haps_sorted_decr_freq.end(),
            [loc_strand] (const pair<const char*,size_t>& a1, const pair<const char*,size_t>& a2) {
                if (a1.second != a2.second)
                    // Decreasing freq.
                    return a1.second > a2.second;
                // Alphabetic (with strand correction).
                if (loc_strand == strand_plus) {
                    return strcmp(a1.first, a2.first) < 0;
                } else {
                    for (const char *s1 = a1.first + strlen(a1.first), *s2 = a2.first + strlen(a2.first);
                            s1 > a1.first && s2 > a2.first;) {
                        --s1; --s2;
                        if (*s1 != *s2)
                            return *s1 > *s2;
                    }
                    return false;
                }
            });

    // Record the indexes.
    for (size_t i=0; i<haps_sorted_decr_freq.size(); i++)
        hap_indexes_map[haps_sorted_decr_freq[i].first] = i;
}

bool Export::is_hap_export()
{
    static const set<std::type_index> hap_exports = {
        typeid(GenePopHapsExport),
        typeid(VcfHapsExport),
        typeid(FineRADStructureExport)
    };
    return hap_exports.count(typeid(*this));
}

int
Export::transpose(ifstream &ifh, vector<string> &transposed)
{
    vector<string> fields;
    string buf;
    size_t len;
    char   line[max_len];

    transposed.clear();

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        ifh.getline(line, max_len);

        //
        // Check if we read a full line.
        //
        while (ifh.fail() && !ifh.eof()) {
            buf += line;
            ifh.clear();
            ifh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        if (!ifh.good() || buf.length() == 0)
            return 0;

        //
        // Break the line up by tabs.
        //
        parse_tsv(buf.c_str(), fields);

        if (transposed.size() == 0) {
            for (uint i = 0; i < fields.size(); i++)
                transposed.push_back(fields[i]);

        } else {
            assert(fields.size() == transposed.size());
            for (uint i = 0; i < fields.size(); i++)
                transposed.at(i) += "\t" + fields[i];
        }

    } while (!ifh.eof());

    return 1;
}

int
SumstatsExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".sumstats.tsv";
    this->_mpopi   = mpopi;
    this->_pop_cnt = this->_mpopi->pops().size();

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cout << "Population-level summary statistics will be written to '" << this->_path << "'\n";

    return 0;
}

int
SumstatsExport::write_header()
{
    //
    // Write the population members.
    //
    for (auto& pop : this->_mpopi->pops()) {
        this->_fh << "# " << pop.name << "\t";
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_fh << this->_mpopi->samples()[i].name;
            if (i < pop.last_sample)
                this->_fh << ",";
        }
        this->_fh << "\n";
    }

    this->_fh << "# Locus ID" << "\t"
              << "Chr"      << "\t"
              << "BP"       << "\t"
              << "Col"      << "\t"
              << "Pop ID"   << "\t"
              << "P Nuc"    << "\t"
              << "Q Nuc"    << "\t"
              << "N"        << "\t"
              << "P"        << "\t"
              << "Obs Het"  << "\t"
              << "Obs Hom"  << "\t"
              << "Exp Het"  << "\t"
              << "Exp Hom"  << "\t"
              << "Pi"       << "\t"
              << "Smoothed Pi"  << "\t"
              << "Smoothed Pi P-value"  << "\t"
              << "Fis"          << "\t"
              << "Smoothed Fis" << "\t"
              << "Smoothed Fis P-value" << "\t"
              << "HWE P-value"  << "\t"
              << "Private"      << "\n";
    return 0;
}

int
SumstatsExport::write_batch(const vector<LocBin *> &loci)
{
    CSLocus        *cloc;
    const LocSum   *s;
    const LocTally *t;
    double          p_freq;

    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;
        t    = loci[i]->s->meta_pop();

        for (uint pos = 0; pos < cloc->len; pos++) {

            if (!t->nucs[pos].fixed) {
                for (uint pop = 0; pop < this->_pop_cnt; pop++) {

                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0)
                        continue;

                    this->_fh
                       << cloc->id << "\t"
                       << cloc->loc.chr() << "\t"
                       << cloc->sort_bp(pos) + 1 << "\t"
                       << pos << "\t"
                       << this->_mpopi->pops()[pop].name << "\t";

                    //
                    // Output the p and q alleles in the same order in each population.
                    //
                    if (t->nucs[pos].p_allele == s->nucs[pos].p_nuc) {
                        if (s->nucs[pos].q_nuc == 0)
                            this->_fh << s->nucs[pos].p_nuc << "\t" << "-";
                        else
                            this->_fh << s->nucs[pos].p_nuc << "\t" << s->nucs[pos].q_nuc;
                        p_freq = s->nucs[pos].p;

                    } else {
                        if (s->nucs[pos].q_nuc == 0)
                            this->_fh << "-\t" << s->nucs[pos].p_nuc;
                        else
                            this->_fh << s->nucs[pos].q_nuc << "\t" << s->nucs[pos].p_nuc;
                        p_freq = 1 - s->nucs[pos].p;
                    }

                    this->_fh << "\t" << (int) s->nucs[pos].num_indv << "\t"
                              << std::setprecision(8)      << p_freq << "\t"
                              << std::setprecision(fieldw) << s->nucs[pos].obs_het << "\t"
                              << s->nucs[pos].obs_hom      << "\t"
                              << s->nucs[pos].exp_het      << "\t"
                              << s->nucs[pos].exp_hom      << "\t"
                              << s->nucs[pos].stat[0]      << "\t" // Pi
                              << s->nucs[pos].smoothed[0]  << "\t"  // Smoothed Pi
                              << s->nucs[pos].bs[0]        << "\t"  // Pi bootstrapped p-value
                              << (s->nucs[pos].stat[1] == -7.0 ? 0.0 : s->nucs[pos].stat[1]) << "\t"  // Fis
                              << s->nucs[pos].smoothed[1]  << "\t"  // Smoothed Fis
                              << s->nucs[pos].bs[1]        << "\t"  // Fis bootstrapped p-value.
                              << s->nucs[pos].stat[2]      << "\t"; // HWE p-value.
                    (t->nucs[pos].priv_allele == (int) pop) ? this->_fh << "1\n" : this->_fh << "0\n";
                }
            }
        }
    }

    return 0;
}

int
HapstatsExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".hapstats.tsv";
    this->_mpopi   = mpopi;
    this->_pop_cnt = this->_mpopi->pops().size();

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cout << "Population-level haplotype summary statistics will be written to '" << this->_path << "'\n";

    return 0;
}

int
HapstatsExport::write_header()
{
    //
    // Write the population members.
    //
    for (auto& pop : this->_mpopi->pops()) {
        this->_fh << "# " << pop.name << "\t";

        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_fh << this->_mpopi->samples()[i].name;
            if (i < pop.last_sample)
                this->_fh << ",";
        }
        this->_fh << "\n";
    }

    this->_fh
        << "# Locus ID"     << "\t"
        << "Chr"            << "\t"
        << "BP"             << "\t"
        << "Pop ID"         << "\t"
        << "N"              << "\t"
        << "Haplotype Cnt"  << "\t"
        << "Gene Diversity" << "\t"
        << "Smoothed Gene Diversity"      << "\t"
        << "Smoothed Gene Diversity P-value"      << "\t"
        << "Haplotype Diversity"          << "\t"
        << "Smoothed Haplotype Diversity" << "\t"
        << "Smoothed Haplotype Diversity P-value" << "\t"
        << "HWE P-value"    << "\t"
        << "HWE P-value SE" << "\t"
        << "Haplotypes"     << "\n";

    return 0;
}

int
HapstatsExport::write_batch(const vector<LocBin *> &loci)
{
    const LocStat *l;

    for (uint i = 0; i < loci.size(); i++) {
        const vector<Pop> &pops = this->_mpopi->pops();

        for (uint pop = 0; pop < this->_pop_cnt; pop++) {

            l = loci[i]->s->hapstats_per_pop(pop);

            if (l == NULL)
                continue;

            this->_fh
                << loci[i]->cloc->id        << "\t"
                << loci[i]->cloc->loc.chr() << "\t"
                << l->bp + 1        << "\t"
                << pops[pop].name   << "\t"
                << (int) l->alleles << "\t"
                << l->hap_cnt       << "\t"
                << l->stat[0]       << "\t"
                << l->smoothed[0]   << "\t"
                << l->bs[0]         << "\t"
                << l->stat[1]       << "\t"
                << l->smoothed[1]   << "\t"
                << l->bs[1]         << "\t"
                << l->stat[2]       << "\t"
                << l->stat[3]       << "\t"
                << l->hap_str       << "\n";
        }
    }

    return 0;
}

SnpDivergenceExport::SnpDivergenceExport(ofstream &log_fh) : _mpopi(NULL), _order(NULL)
{
    if (ordered_export)
        this->_order = new OPopPair<PopPair>(log_fh);
}

int
SnpDivergenceExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Open an SNP Divergence output file for each pair of populations.
    //
    string    path;
    ofstream *fh;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        const Pop& pop_1p = this->_mpopi->pops()[pop_1];

        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            const Pop& pop_2p = this->_mpopi->pops()[pop_2];

            path = out_path + out_prefix + ".fst_" + pop_1p.name + "-" + pop_2p.name + ".tsv";
            fh = new ofstream(path.c_str(), ofstream::out);
            check_open(*fh, path);
            fh->precision(fieldw);
            fh->setf(std::ios::fixed);

            this->_fhs.push_back(fh);
        }
    }

    return 0;
}

int
SnpDivergenceExport::write_header()
{
    //
    // Write a header into each SNP Divergence output file (one for each pair of populations).
    //
    ofstream *fh;

    uint i = 0;
    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            fh = this->_fhs[i];

            *fh << "# Locus ID" << "\t"
                << "Pop 1 ID"   << "\t"
                << "Pop 2 ID"   << "\t"
                << "Chr"        << "\t"
                << "BP"         << "\t"
                << "Column"     << "\t"
                << "Overall Pi" << "\t"
                << "AMOVA Fst"  << "\t"
                << "Fisher's P" << "\t"
                << "Odds Ratio" << "\t"
                << "CI Low"     << "\t"
                << "CI High"    << "\t"
                << "LOD"        << "\t"
                << "Corrected AMOVA Fst"        << "\t"
                << "Smoothed AMOVA Fst"         << "\t"
                << "Smoothed AMOVA Fst P-value" << "\t"
                << "Window SNP Count";

            //
            // If requested, log Fst component calculations to a file.
            //
            if (log_fst_comp) {
                *fh << "\t"
                    << "n_1" << "\t"
                    << "n_2" << "\t"
                    << "tot_alleles" << "\t"
                    << "p_1" << "\t"
                    << "q_1" << "\t"
                    << "p_2" << "\t"
                    << "q_2" << "\t"
                    << "pi_1" << "\t"
                    << "pi_2" << "\t"
                    << "pi_all" << "\t"
                    << "bcoeff_1" << "\t"
                    << "bcoeff_2" << "\t"
                    << "binomial_fst" << "\t"
                    << "p_1_freq" << "\t"
                    << "q_1_freq" << "\t"
                    << "p_2_freq" << "\t"
                    << "q_2_freq" << "\t"
                    << "p_avg_cor" << "\t"
                    << "n_avg_cor" << "\t"
                    << "amova_fst" << "\n";
            } else {
                *fh << "\n";
            }

            i++;
        }
    }
    return 0;
}

int
SnpDivergenceExport::write_batch_pairwise(const vector<LocBin *> &loci, const vector<vector<PopPair **>> &div)
{
    vector<const PopPair *> sites;

    for (uint i = 0; i < div.size(); i++) {

        assert(div[i].size() == loci.size());
        sites.clear();

        if (ordered_export) {
            this->_order->order(sites, loci, div[i]);

            string chr = loci[0]->cloc->loc.chr();

            for (uint j = 0; j < sites.size(); j++)
                if (sites[j] != NULL)
                    this->write_site(this->_fhs[i], sites[j], chr);

        } else {
            PopPair **pp;
            size_t    cloc_len;
            LocBin   *loc;

            for (uint j = 0; j < div[i].size(); j++) {

                loc = loci[j];
                cloc_len = loc->cloc->len;
                pp       = div[i][j];

                for (uint pos = 0; pos < cloc_len; pos++)
                    if (pp[pos] != NULL)
                        this->write_site(this->_fhs[i], pp[pos], loc->cloc->loc.chr());
            }
        }
    }

    return 0;
}

inline int
SnpDivergenceExport::write_site(ofstream *fh, const PopPair *pp, string chr)
{
    *fh << pp->loc_id      << "\t"
        << this->_mpopi->pops()[pp->pop_1].name << "\t"
        << this->_mpopi->pops()[pp->pop_2].name << "\t"
        << chr             << "\t"
        << pp->bp + 1      << "\t"
        << pp->col         << "\t"
        << pp->pi          << "\t"
        << pp->amova_fst   << "\t"
        << std::setprecision(9)      << pp->fet_p  << "\t"
        << std::setprecision(fieldw) << pp->fet_or << "\t"
        << pp->ci_low      << "\t"
        << pp->ci_high     << "\t"
        << pp->lod         << "\t"
        << pp->stat[1]     << "\t"
        << pp->smoothed[1] << "\t"
        << pp->bs[1]       << "\t"
        << pp->snp_cnt;

    if (log_fst_comp) {
        *fh << "\t"
            << pp->comp[0]   << "\t"
            << pp->comp[1]   << "\t"
            << pp->comp[2]   << "\t"
            << pp->comp[3]   << "\t"
            << pp->comp[4]   << "\t"
            << pp->comp[5]   << "\t"
            << pp->comp[6]   << "\t"
            << pp->comp[7]   << "\t"
            << pp->comp[8]   << "\t"
            << pp->comp[9]   << "\t"
            << pp->comp[10]  << "\t"
            << pp->comp[11]  << "\t"
            << pp->fst       << "\t"
            << pp->comp[12]  << "\t"
            << pp->comp[13]  << "\t"
            << pp->comp[14]  << "\t"
            << pp->comp[15]  << "\t"
            << pp->comp[16]  << "\t"
            << pp->comp[17]  << "\t"
            << pp->amova_fst << "\n";
    } else {
        *fh << "\n";
    }

    return 0;
}

int
HapDivergenceExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Open an Haplotype Divergence output file for each pair of populations.
    //
    string    path;
    ofstream *fh;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        const Pop& pop_1p = this->_mpopi->pops()[pop_1];

        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            const Pop& pop_2p = this->_mpopi->pops()[pop_2];

            path = out_path + out_prefix + ".phistats_" + pop_1p.name + "-" + pop_2p.name + ".tsv";
            fh = new ofstream(path.c_str(), ofstream::out);
            check_open(*fh, path);
            fh->precision(fieldw);
            fh->setf(std::ios::fixed);

            this->_fhs.push_back(fh);
        }
    }

    //
    // Open an single Haplotype Divergence output file for the full meta population.
    //
    path = out_path + out_prefix + ".phistats.tsv";
    fh = new ofstream(path.c_str(), ofstream::out);
    check_open(*fh, path);
    fh->precision(fieldw);
    fh->setf(std::ios::fixed);

    this->_metapop_fh = fh;

    return 0;
}

int
HapDivergenceExport::write_header()
{
    //
    // Write a header into each haplotype divergence output file (one for each pair of populations).
    //
    ofstream *fh;

    uint i = 0;
    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            vector<int> subpop_ids;
            subpop_ids.push_back(pop_1);
            subpop_ids.push_back(pop_2);
            fh = this->_fhs[i];

            //
            // Write the population members.
            //
            for (int k : subpop_ids) {
                const Pop& pop_k = this->_mpopi->pops()[k]; // This is [pop_i], then [pop_j].
                *fh << "# Population " << pop_k.name << "\t";
                for (size_t n = pop_k.first_sample; n <= pop_k.last_sample; n++) {
                    *fh << this->_mpopi->samples()[n].name;
                    if (n < pop_k.last_sample)
                        *fh << ",";
                }
                *fh << "\n";
            }

            *fh << "# Locus ID" << "\t"
                << "Pop 1 ID"   << "\t"
                << "Pop 2 ID"   << "\t"
                << "Chr"        << "\t"
                << "BP"         << "\t";
            if (log_fst_comp)
                *fh << "SSD(WP)"     << "\t"
                    << "SSD(AP/WG)"  << "\t"
                    << "SSD(AG)"     << "\t"
                    << "SSD(TOTAL)"  << "\t"
                    << "MSD(WP)"     << "\t"
                    << "MSD(AP/WG)"  << "\t"
                    << "MSD(AG)"     << "\t"
                    << "MSD(TOTAL)"  << "\t"
                    << "n"           << "\t"
                    << "n'"          << "\t"
                    << "n''"         << "\t"
                    << "Sigma2_a"    << "\t"
                    << "Sigma2_b"    << "\t"
                    << "Sigma2_c"    << "\t"
                    << "Sigma_Total" << "\t";
            *fh << "phi_st"          << "\t"
                << "Smoothed Phi_st" << "\t"
                << "Smoothed Phi_st P-value" << "\t"
                << "Fst'"            << "\t"
                << "Smoothed Fst'"   << "\t"
                << "Smoothed Fst' P-value"   << "\t"
                << "D_est"          << "\t"
                << "Smoothed D_est" << "\t"
                << "Smoothed D_est P-value" << "\t"
                << "Dxy"            << "\t"
                << "Smoothed Dxy"   << "\t"
                << "Smoothed Dxy P-value"   << "\n";

            i++;
        }
    }

    //
    // Write a header into the meta population haplotype divergence output file.
    //
    fh = this->_metapop_fh;

    //
    // Write the population members.
    //
    for (auto& pop : this->_mpopi->pops()) {
        *fh << "# Population " << pop.name << "\t";
        for (size_t k = pop.first_sample; k <= pop.last_sample; k++) {
            *fh << this->_mpopi->samples()[k].name;
            if (k < pop.last_sample)
                *fh << ",";
        }
        *fh << "\n";
    }

    //
    // Write the group members.
    //
    for (auto& group : this->_mpopi->groups()) {
        *fh << "# Group " << group.name << "\t";
        for (size_t i_pop : group.pops) {
            *fh << this->_mpopi->pops()[i_pop].name;
            if (i_pop != group.pops.back())
                *fh << ",";
        }
        *fh << "\n";
    }

    *fh << "# Locus ID" << "\t"
        << "Chr"        << "\t"
        << "BP"         << "\t";
    if (log_fst_comp)
        *fh << "SSD(WP)"     << "\t"
            << "SSD(AP/WG)"  << "\t"
            << "SSD(AG)"     << "\t"
            << "SSD(TOTAL)"  << "\t"
            << "MSD(WP)"     << "\t"
            << "MSD(AP/WG)"  << "\t"
            << "MSD(AG)"     << "\t"
            << "MSD(TOTAL)"  << "\t"
            << "n"           << "\t"
            << "n'"          << "\t"
            << "n''"         << "\t"
            << "Sigma2_a"    << "\t"
            << "Sigma2_b"    << "\t"
            << "Sigma2_c"    << "\t"
            << "Sigma_Total" << "\t";
    *fh << "phi_st"          << "\t"
        << "Smoothed Phi_st" << "\t"
        << "Smoothed Phi_st P-value" << "\t"
        << "Fst'"            << "\t"
        << "Smoothed Fst'"   << "\t"
        << "Smoothed Fst' P-value"   << "\t"
        << "D_est"          << "\t"
        << "Smoothed D_est" << "\t"
        << "Smoothed D_est P-value" << "\n";

    return 0;
}

int
HapDivergenceExport::write_batch_pairwise(const vector<LocBin *> &loci,
                                          const vector<vector<HapStat *>> &div,
                                          const vector<HapStat *> &metapop_div)
{
    ofstream *fh;
    HapStat  *h;
    LocBin   *loc;

    for (uint i = 0; i < div.size(); i++) {

        fh = this->_fhs[i];

        assert(div[i].size() == loci.size());

        for (uint j = 0; j < div[i].size(); j++) {
            loc = loci[j];
            h   = div[i][j];

            if (h == NULL) continue;

            *fh << h->loc_id    << "\t"
                << this->_mpopi->pops()[h->pop_1].name << "\t"
                << this->_mpopi->pops()[h->pop_2].name << "\t"
                << loc->cloc->loc.chr() << "\t"
                << h->bp +1     << "\t";
            if (log_fst_comp)
                *fh << h->comp[0]  << "\t"
                    << h->comp[1]  << "\t"
                    << h->comp[2]  << "\t"
                    << h->comp[3]  << "\t"
                    << h->comp[4]  << "\t"
                    << h->comp[5]  << "\t"
                    << h->comp[6]  << "\t"
                    << h->comp[7]  << "\t"
                    << h->comp[8]  << "\t"
                    << h->comp[9]  << "\t"
                    << h->comp[10] << "\t"
                    << h->comp[11] << "\t"
                    << h->comp[12] << "\t"
                    << h->comp[13] << "\t"
                    << h->comp[14] << "\t";
            *fh << h->stat[0]     << "\t"
                << h->smoothed[0] << "\t"
                << h->bs[0]       << "\t"
                << h->stat[3]     << "\t"
                << h->smoothed[3] << "\t"
                << h->bs[3]       << "\t"
                << h->stat[4]     << "\t"
                << h->smoothed[4] << "\t"
                << h->bs[4]       << "\t"
                << h->stat[5]     << "\t"
                << h->smoothed[5] << "\t"
                << h->bs[5]       << "\n";
        }
    }

    fh = this->_metapop_fh;

    assert(metapop_div.size() == loci.size());

    for (uint j = 0; j < metapop_div.size(); j++) {
        loc = loci[j];
        h   = metapop_div[j];

        if (h == NULL) continue;

        *fh << h->loc_id    << "\t"
            << loc->cloc->loc.chr() << "\t"
            << h->bp +1     << "\t";
        if (log_fst_comp)
            *fh << h->comp[0]  << "\t"
                << h->comp[1]  << "\t"
                << h->comp[2]  << "\t"
                << h->comp[3]  << "\t"
                << h->comp[4]  << "\t"
                << h->comp[5]  << "\t"
                << h->comp[6]  << "\t"
                << h->comp[7]  << "\t"
                << h->comp[8]  << "\t"
                << h->comp[9]  << "\t"
                << h->comp[10] << "\t"
                << h->comp[11] << "\t"
                << h->comp[12] << "\t"
                << h->comp[13] << "\t"
                << h->comp[14] << "\t";
        *fh << h->stat[0]     << "\t"
            << h->smoothed[0] << "\t"
            << h->bs[0]       << "\t"
            << h->stat[3]     << "\t"
            << h->smoothed[3] << "\t"
            << h->bs[3]       << "\t"
            << h->stat[4]     << "\t"
            << h->smoothed[4] << "\t"
            << h->bs[4]       << "\n";
    }

    return 0;
}

int
RawHaplotypesExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".haplotypes.tsv";
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    cout << "Raw haplotypes will be written to '" << this->_path << "'\n";

    return 0;
}

int
RawHaplotypesExport::write_header()
{
    this->_fh << "# Catalog Locus ID" << "\t" << "Cnt";

    for (size_t i : this->_mpopi->sample_indexes_orig_order())
        this->_fh << "\t" << this->_mpopi->samples()[i].name;
    this->_fh << "\n";

    return 0;
}

int
RawHaplotypesExport::write_batch(const vector<LocBin *> &loci)
{
    CSLocus *cloc;

    //
    // Output each locus.
    //
    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;

        stringstream id;
        cloc->annotation.length() > 0 ?
            id << cloc->id << "|" << cloc->annotation : id << cloc->id;

        this->_fh << id.str();

        this->_fh << "\t" << cloc->cnt;

        Datum **d = loci[i]->d;
        string  obshap;

        for (size_t i : this->_mpopi->sample_indexes_orig_order()) {
            this->_fh << "\t";

            if (d[i] == NULL)
                this->_fh << "-";
            else {
                obshap = "";
                for (uint j = 0; j < d[i]->obshap.size(); j++)
                    obshap += string(d[i]->obshap[j]) + "/";
                obshap = obshap.substr(0, obshap.length()-1);
                this->_fh << obshap;
            }
        }

        this->_fh << "\n";
    }

    return 0;
}

int
MarkersExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".markers.tsv";
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cout << "Genotyping markers will be written to '" << this->_path << "'\n";

    return 0;
}

int
MarkersExport::write_header()
{
    this->_fh << "# LocusID\t"
              << "Marker\t"
              << "Genotypes\t"
              << "MappableGenotypes\t"
              << "Segregation Distortion P-Value\n";
    return 0;
}

int
MarkersExport::write_batch(const vector<LocBin *> &loci)
{
    stringstream gtype_map;

    for (uint i = 0; i < loci.size(); i++) {
        this->_fh
            << loci[i]->cloc->id     << "\t"
            << loci[i]->cloc->marker << "\t" 
            << loci[i]->cloc->cnt    << "\t"
            << loci[i]->cloc->gcnt   << "\t"
            << loci[i]->cloc->chisq 
            << "\n";
    }

    return 0;
}

int
FastaLociExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".loci.fa";
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    cout << "FASTA consensus sequences for each locus in the metapopulation  will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaLociExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "# Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaLociExport::write_batch(const vector<LocBin *> &loci)
{
    for (uint i = 0; i < loci.size(); i++) {
        LocBin* loc = loci[i];

        this->_fh << ">CLocus_" << loc->cloc->id;
        if (strcmp(loc->cloc->loc.chr(), "un") != 0)
            this->_fh << " [" << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-") << "]";
        this->_fh << '\n' << loc->cloc->con << '\n';

#ifdef DEBUG
        bool no_samples = true;
        Datum** d = loc->d;
        for (uint j = 0; j < this->_mpopi->samples().size(); j++)
            if (d[j] != NULL)
                no_samples = false;
        assert(!no_samples);
#endif
    }

    return 0;
}

int
FastaRawExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "# Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaRawExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_path = out_path + out_prefix + ".samples-raw.fa";
    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    cout << "Raw FASTA consensus sequences for each sample will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaRawExport::write_batch(const vector<LocBin *> &loci)
{
    //
    // Write a FASTA file containing each allele from each locus from
    // each sample in the population.
    //
    LocBin *loc;
    Datum **d;
    char   *seq;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i];
        d   = loc->d;
        seq = new char[loc->cloc->len + 1];
        strcpy(seq, loc->cloc->con);

        for (uint j : this->_mpopi->sample_indexes_orig_order()) {
            if (d[j] == NULL)
                continue;

            for (uint k = 0; k < d[j]->obshap.size(); k++) {

                for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                    uint col = loc->cloc->snps[i]->col;
                    seq[col] = col < loc->cloc->len ? d[j]->obshap[k][i] : loc->cloc->con[col];
                }

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << k
                          << " ["       << this->_mpopi->samples()[j].name;

                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";
            }
        }
        delete [] seq;
    }

    return 0;
}

int
FastaSamplesExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "# Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaSamplesExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_path = out_path + out_prefix + ".samples.fa";
    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    cout << "FASTA consensus sequences for each sample will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaSamplesExport::write_batch(const vector<LocBin *> &loci)
{
    LocBin *loc;
    Datum **d;
    char   *seq;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i];
        d   = loc->d;
        seq = new char[loc->cloc->len + 1];
        strcpy(seq, loc->cloc->con);

        for (uint j : this->_mpopi->sample_indexes_orig_order()) {
            if (d[j] == NULL)
                continue;
            if (d[j]->obshap.size() > 2)
                continue;

            if (d[j]->obshap.size() == 1) {

                for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                    uint col = loc->cloc->snps[i]->col;
                    seq[col] = col < loc->cloc->len ? d[j]->obshap[0][i] : loc->cloc->con[col];
                }

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << 0
                          << " ["       << this->_mpopi->samples()[j].name;
                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << 1
                          << " ["       << this->_mpopi->samples()[j].name;
                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";

            } else {
                for (uint k = 0; k < d[j]->obshap.size(); k++) {
                    for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                        uint col = loc->cloc->snps[i]->col;
                        seq[col] = col < loc->cloc->len ? d[j]->obshap[k][i] : loc->cloc->con[col];
                    }

                    this->_fh << ">CLocus_" << loc->cloc->id
                              << "_Sample_" <<  this->_mpopi->samples()[j].id
                              << "_Locus_"  << d[j]->id
                              << "_Allele_" << k
                              << " ["       <<  this->_mpopi->samples()[j].name;
                    if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                        this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                    this->_fh << "]\n"
                              << seq << "\n";
                }
            }
        }

        delete [] seq;
    }

    return 0;
}

int
OrderableExport::write_batch(const vector<LocBin*> &loci)
{
    if (ordered_export) {
        //
        // We need to order the SNPs to take into account overlapping loci.
        //
        vector<const NucTally *> sites;
        OLocTally<NucTally> ord;
        ord.order(sites, loci);

        map<size_t, size_t> key;
        map<size_t, size_t>::iterator key_it = key.begin();

        for (size_t k = 0; k < loci.size(); k++)
            key_it = key.insert(key_it, pair<size_t, size_t>(loci[k]->cloc->id, k));

        for (uint pos = 0; pos < sites.size(); pos++) {
            int loc_id = sites[pos]->loc_id;

            const LocBin* loc = loci[key[loc_id]];
            size_t        col = sites[pos]->col;
            size_t  snp_index = loc->cloc->snp_index(col);

            this->write_site(loc->cloc, loc->s, loc->d, col, snp_index);
        }

    } else {
        for (uint k = 0; k < loci.size(); k++) {
            const LocBin*    loc = loci[k];
            const CSLocus*  cloc = loc->cloc;
            Datum const*const* d = loc->d;
            const LocTally*    t = loc->s->meta_pop();

            for (uint snp_index = 0; snp_index < cloc->snps.size(); snp_index++) {
                uint col = cloc->snps[snp_index]->col;

                if (t->nucs[col].allele_cnt != 2)
                    continue;

                this->write_site(cloc, loc->s, d, col, snp_index);
            }
        }
    }

    return 0;
}

int
GenePopExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Write a GenePop file as defined here: http://kimura.univ-montp2.fr/~rousset/Genepop.htm
    //
    this->_path = out_path + out_prefix + ".snps.genepop";

    //
    // Open a temporary file.
    //
    this->_tmpfh.open(this->tmp_path());
    check_open(this->_tmpfh, this->tmp_path());

    cout << "Polymorphic sites in GenePop format will be written to '" << this->_path << "'\n";

    return 0;
}

int
GenePopExport::write_header()
{
    //
    // N.B. For the GENEPOP output, we iterate over samples in the sorted order
    // rather than in the original one, because samples MUST be grouped by
    // population.
    //
    const vector<size_t>& orig_order = this->_mpopi->sample_indexes_orig_order();
    if (!std::is_sorted(orig_order.begin(), orig_order.end()))
        cerr << "Warning: Genepop: The order in which samples appear was modified"
            " (as the input population map is not sorted).";
    for (const Pop& pop : this->_mpopi->pops())
        for (size_t j = pop.first_sample; j <= pop.last_sample; j++)
            this->_tmpfh << "\t" << this->_mpopi->samples()[j].name;
    this->_tmpfh << "\n";

    return 0;
}

int
GenePopExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    map<char, string> nuc_map;
    nuc_map['A'] = "01";
    nuc_map['C'] = "02";
    nuc_map['G'] = "03";
    nuc_map['T'] = "04";

    const LocSum *s;
    char p_allele, q_allele;

    this->_tmpfh << loc->id << "_" << col;

    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        s = lps->per_pop(p);

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {

            if (s->nucs[col].incompatible_site ||
                s->nucs[col].filtered_site) {
                //
                // This site contains more than two alleles in this population or was filtered
                // due to a minor allele frequency that is too low.
                //
                this->_tmpfh << "\t0000";
            } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                //
                // Data does not exist.
                //
                this->_tmpfh << "\t0000";
            } else if (d[j]->model[col] == 'U') {
                //
                // Data exists, but the model call was uncertain.
                //
                this->_tmpfh << "\t0000";
            } else {
                //
                // Tally up the nucleotide calls.
                //
                tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

                if (p_allele == 0 && q_allele == 0) {
                    // More than two potential alleles.
                    this->_tmpfh << "\t0000";
                } else if (p_allele == 0) {
                    this->_tmpfh << "\t" << nuc_map[q_allele] << nuc_map[q_allele];

                } else if (q_allele == 0) {
                    this->_tmpfh << "\t" << nuc_map[p_allele] << nuc_map[p_allele];

                } else {
                    this->_tmpfh << "\t" << nuc_map[p_allele] << nuc_map[q_allele];
                }
            }
        }
    }
    this->_tmpfh << "\n";
    return 0;
}

int
GenePopExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh.open(this->_path.c_str());
    check_open(this->_fh, this->_path);

    //
    // Output the header line.
    //
    this->_fh << "# Stacks v" << VERSION << "; GenePop v4.1.3; " << date << "\n";

    ifstream intmpfh (this->tmp_path());
    check_open(intmpfh, this->tmp_path());

    vector<string> transposed_lines;

    Export::transpose(intmpfh, transposed_lines);

    assert(transposed_lines.size() == this->_mpopi->samples().size() + 1);

    size_t line_cnt = 0;
    size_t pos;
    //
    // The first line has a list of locus IDs, convert these to comma-separated.
    //
    for (uint i = 1; i < transposed_lines[line_cnt].size(); i++)
        if (transposed_lines[line_cnt][i] == '\t')
            transposed_lines[line_cnt][i] = ',';

    this->_fh << transposed_lines[line_cnt].substr(1) << "\n";
    line_cnt++;

    for (const Pop& pop : this->_mpopi->pops()) {
        this->_fh << "pop\n";
        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
            pos = transposed_lines[line_cnt].find_first_of('\t');
            transposed_lines[line_cnt][pos] = ',';
            this->_fh
                << transposed_lines[line_cnt].substr(0, pos + 1) << "\t"
                << transposed_lines[line_cnt].substr(pos + 1) << "\n";
            line_cnt++;
        }
    }

    return 1;
}

int
GenePopHapsExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;
    this->_path = out_path + out_prefix + ".haps.genepop";
    this->_fh.open(this->_path);
    check_open(this->_fh, this->_path);
    this->_tmpfh.open(this->tmp_path());
    check_open(this->_tmpfh, this->tmp_path());
    cout << "Polymorphic loci in GenePop format will be written to '" << this->_path << "'\n";
    return 0;
}

int
GenePopHapsExport::write_header()
{
    for (const Pop& pop : this->_mpopi->pops()) // Per pop; c.f. GenePopExport.
        for (size_t s = pop.first_sample; s <= pop.last_sample; s++)
            this->_tmpfh << "\t" << this->_mpopi->samples()[s].name << ',';
    this->_tmpfh << "\n";
    this->_tmpfh << std::setfill('0');
    return 0;
}

int
GenePopHapsExport::write_batch(const vector<LocBin*>& loci)
{
    for (const LocBin* locbin : loci) {
        const CSLocus* cloc = locbin->cloc;
        Datum const*const* d = locbin->d;

        if (cloc->snps.empty())
            // Monomorphic locus.
            continue;

        // Tally and sort haplotypes.
        vector<pair<const char*, size_t>> sorted_haps;
        map<const char*, size_t, LessCStrs> hap_indexes;
        tally_complete_haplotypes(d, this->_mpopi->samples().size(), cloc->loc.strand,
                                  sorted_haps, hap_indexes);
        if (sorted_haps.size() < 2)
            continue;

        this->_tmpfh << cloc->id;
        for (const Pop& pop : this->_mpopi->pops()) {
            for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                this->_tmpfh << "\t";
                if (d[j] == NULL || strchr(d[j]->obshap[0], 'N') != NULL) {
                    this->_tmpfh << (_n_digits == 2 ? "0000" : "000000");
                } else {
                    size_t i0 = hap_indexes.at(d[j]->obshap[0]);
                    size_t i1 = hap_indexes.at(d[j]->obshap[1]);
                    if (max(i0, i1) + 1 > (_n_digits == 2 ? 99 : 999)) {
                        this->_tmpfh << (_n_digits == 2 ? "0000" : "000000");
                        continue;
                    }
                    this->_tmpfh << std::setw(_n_digits) << min(i0, i1) + 1
                                 << std::setw(_n_digits) << max(i0, i1) + 1;
                }
            }
        }
        this->_tmpfh << "\n";
    }
    return 0;
}

int
GenePopHapsExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header line.
    //
    this->_fh << "Stacks v" << VERSION << "; GenePop v4.1.3; " << date << "\n";

    //
    // Read and transpose the temporary file.
    //
    ifstream intmpfh (this->tmp_path());
    check_open(intmpfh, this->tmp_path());
    vector<string> transposed_lines;
    Export::transpose(intmpfh, transposed_lines);
    if (transposed_lines.empty()) {
        cerr << "Error: Temporary file '" << this->tmp_path()
             << "' is corrupt (no data).\n";
        throw exception();
    } else if (transposed_lines.size() != this->_mpopi->samples().size() + 1) {
        cerr << "Error: Temporary file '" << this->tmp_path()
             << "' is corrupt (wrong number of columns).\n";
        throw exception();
    }

    //
    // The first line has a list of locus IDs, convert these to comma-separated.
    //
    vector<string>::iterator line = transposed_lines.begin();
    for (auto ch = ++line->begin(); ch != line->end(); ++ch) {
        if (*ch == '\t')
            this->_fh << ',';
        else
            this->_fh << *ch;
    }
    this->_fh << '\n';
    line++;

    //
    // Output every sample.
    //
    for (const Pop& pop : this->_mpopi->pops()) {
        this->_fh << "pop\n";
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++, line++) {
            assert(line != transposed_lines.end());
            this->_fh << *line << '\n';
        }
    }
    assert(line == transposed_lines.end());

    return 0;
}

int
StructureExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Write a Structure file as defined here: http://pritch.bsd.uchicago.edu/structure.html
    //
    this->_path = out_path + out_prefix + ".structure";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + ".structure.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "Polymorphic sites in Structure format will be written to '" << this->_path << "'\n";

    return 0;
}

int
StructureExport::write_header()
{
    for (size_t j : this->_mpopi->sample_indexes_orig_order())
        this->_tmpfh << '\t' << this->_mpopi->samples()[j].name << '\t' << this->_mpopi->samples()[j].name;
    this->_tmpfh << '\n';

    for (size_t j : this->_mpopi->sample_indexes_orig_order()) {
        const Sample& sample = this->_mpopi->samples()[j];
        const Pop& pop = this->_mpopi->pops()[sample.pop];
        this->_tmpfh << '\t' << pop.name << '\t' << pop.name;
    }
    this->_tmpfh << "\n";

    return 0;
}

int
StructureExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    map<char, string> nuc_map;
    nuc_map['A'] = "1";
    nuc_map['C'] = "2";
    nuc_map['G'] = "3";
    nuc_map['T'] = "4";

    const LocSum *s;
    char p_allele, q_allele;

    this->_tmpfh << loc->id << "_" << col;

    for (size_t j : this->_mpopi->sample_indexes_orig_order()) {
        s = lps->per_pop(this->_mpopi->samples()[j].pop);
        if (s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            //
            // This site contains more than two alleles in this population or was filtered
            // due to a minor allele frequency that is too low.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
            //
            // Data does not exist.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j]->model[col] == 'U') {
            //
            // Data exists, but the model call was uncertain.
            //
            this->_tmpfh << "\t" << "0";
        } else {
            //
            // Tally up the nucleotide calls.
            //
            tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

            if (p_allele == 0 && q_allele == 0)
                this->_tmpfh << "\t" << "0";
            else if (p_allele == 0)
                this->_tmpfh << "\t" << nuc_map[q_allele];
            else
                this->_tmpfh << "\t" << nuc_map[p_allele];
        }

        //
        // Output the site for this sample again, now for the q allele
        //
        if (s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            this->_tmpfh << "\t" << "0";
        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
            this->_tmpfh << "\t" << "0";
        } else if (d[j]->model[col] == 'U') {
            this->_tmpfh << "\t" << "0";
        } else {
            tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

            if (p_allele == 0 && q_allele == 0)
                this->_tmpfh << "\t" << "0";
            else if (q_allele == 0)
                this->_tmpfh << "\t" << nuc_map[p_allele];
            else
                this->_tmpfh << "\t" << nuc_map[q_allele];
        }
    }
    this->_tmpfh << "\n";

    return 0;
}

int
StructureExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Output the header line.
    //
    this->_fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n";

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    vector<string> transposed_lines;

    Export::transpose(this->_intmpfh, transposed_lines);

    assert(transposed_lines.size() == (this->_mpopi->samples().size() * 2) + 1);

    for (size_t line_cnt = 0; line_cnt < transposed_lines.size(); line_cnt++)
        this->_fh << transposed_lines[line_cnt] << "\n";

    return 1;
}

void
StructureExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();

    remove(this->_tmp_path.c_str());

    this->_fh.close();
    return;
}

int
FineRADStructureExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;
    this->_path = out_path + out_prefix + ".haps.radpainter";
    cout << "Polymorphic loci in RADpainter/fineRADstructure format will be written to '"
         << this->_path << "'\n";
    this->_fh.open(this->_path.c_str());
    check_open(this->_fh, this->_path);
    return 0;
}

int
FineRADStructureExport::write_header()
{
    bool first = true;
    for (size_t sample : this->_mpopi->sample_indexes_orig_order()) {
        if (first)
            first = false;
        else
            this->_fh << '\t';
        this->_fh << this->_mpopi->samples()[sample].name;
    }
    this->_fh << '\n';
    return 0;
}

int
FineRADStructureExport::write_batch(const vector<LocBin *> &loci)
{
    const CSLocus*     cloc;
    Datum const*const* d;

    for (const LocBin* loc : loci) {
        cloc = loc->cloc;
        d    = loc->d;

        // Monomorphic locus.
        if (cloc->snps.empty())
            continue;

        //
        // Tally and sort haplotypes.
        //
        vector<pair<const char*, size_t>>   sorted_haps;
        map<const char*, size_t, LessCStrs> hap_indexes;
        tally_complete_haplotypes(d, this->_mpopi->samples().size(), cloc->loc.strand,
                                  sorted_haps, hap_indexes);
        if (sorted_haps.size() < 2)
            continue;

        bool first = true;
        for (size_t s : this->_mpopi->sample_indexes_orig_order()) {
            if (first)
                first = false;
            else
                this->_fh << '\t';
            if (d[s] == NULL || strchr(d[s]->obshap[0], 'N') != NULL)
                continue;
            int cmp = strcmp(d[s]->obshap[0], d[s]->obshap[1]);
            if (cmp == 0)
                this->_fh << d[s]->obshap[0];
            else if (cmp < 0)
                this->_fh << d[s]->obshap[0] << '/' << d[s]->obshap[1];
            else
                this->_fh << d[s]->obshap[1] << '/' << d[s]->obshap[0];
        }
        this->_fh << '\n';
    }

    return 0;
}

inline char
iupac_encode(char p_nuc, char q_nuc)
{
    char nuc = '\0';
    //
    // Encode SNPs that are variable within a population using IUPAC notation:
    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
    //
    switch(p_nuc) {
    case 0:
        nuc = 'N';
        break;
    case 'A':
        switch(q_nuc) {
        case 'C':
            nuc = 'M';
            break;
        case 'G':
            nuc = 'R';
            break;
        case 'T':
            nuc = 'W';
            break;
        case 0:
            nuc = 'A';
            break;
        }
        break;
    case 'C':
        switch(q_nuc) {
        case 'A':
            nuc = 'M';
            break;
        case 'G':
            nuc = 'S';
            break;
        case 'T':
            nuc = 'Y';
            break;
        case 0:
            nuc = 'C';
            break;
        }
        break;
    case 'G':
        switch(q_nuc) {
        case 'A':
            nuc = 'R';
            break;
        case 'C':
            nuc = 'S';
            break;
        case 'T':
            nuc = 'K';
            break;
        case 0:
            nuc = 'G';
            break;
        }
        break;
    case 'T':
        switch(q_nuc) {
        case 'A':
            nuc = 'W';
            break;
        case 'C':
            nuc = 'Y';
            break;
        case 'G':
            nuc = 'K';
            break;
        case 0:
            nuc = 'T';
            break;
        }
        break;
    }

    return nuc;
}

int
PhylipExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    string suffix = typeid(*this) == typeid(PhylipVarExport) ? ".var" : ".fixed";

    //
    // We will write those loci to a Phylip file as defined here:
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    this->_path = out_path + out_prefix + suffix + ".phylip";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + suffix + ".phylip.part";
    this->_tmpfh.open(this->_tmp_path);

    //
    // Open a log file.
    //
    this->_log_path = this->_path + ".log";
    this->_logfh.open(this->_log_path);

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_logfh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n"
                 << "# Seq Pos\tLocus ID\tColumn\tPopulation\n";

    if (typeid(*this) == typeid(PhylipVarExport))
        cout << "Polymorphic";
    else
        cout << "Fixed difference";
    cout << " sites in Phylip format will be written to '" << this->_path << "'\n";

    return 0;
}

int
PhylipExport::write_header()
{
    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        this->_tmpfh << pop.name;
        if (p < this->_mpopi->pops().size() - 1)
            this->_tmpfh << "\t";
    }
    this->_tmpfh << "\n";

    return 0;
}

int
PhylipFixedExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    const LocSum   *s;
    const LocTally *t = lps->meta_pop();
    size_t pop_cnt = this->_mpopi->pops().size();

    //
    // We are looking for loci that are fixed within each population, but are
    // variable between one or more populations.
    //
    if (t->nucs[col].fixed == true || t->nucs[col].allele_cnt != 2 || t->nucs[col].pop_cnt < 2)
        return 0;

    bool fixed_within = true;
    for (uint p = 0; p < pop_cnt; p++) {
        s = lps->per_pop(p);

        if (s->nucs[col].num_indv == 0)
            continue;
        if (s->nucs[col].fixed == false) {
            fixed_within = false;
            break;
        }
    }
    if (fixed_within == false)
        return 0;

    this->_logfh << this->_site_index << "\t" << loc->id << "\t" << col << "\t";

    for (uint p = 0; p < pop_cnt; p++) {
        s = lps->per_pop(p);

        if (s->nucs[col].num_indv > 0) {
            this->_tmpfh << s->nucs[col].p_nuc;
            this->_logfh << this->_mpopi->pops()[p].name << ":" << s->nucs[col].p_nuc << ",";
        } else {
            this->_tmpfh << "N";
            this->_logfh << this->_mpopi->pops()[p].name << ":N" << ",";
        }
        if (p < pop_cnt - 1) this->_tmpfh << "\t";
    }
    this->_logfh << "\n";
    this->_site_index++;

    this->_tmpfh << "\n";

    return 0;
}

int
PhylipVarExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    const LocSum   *s;
    const LocTally *t = lps->meta_pop();
    size_t pop_cnt = this->_mpopi->pops().size();

    //
    // Encode SNPs that are variable within a population as well, using IUPAC notation:
    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
    //
    if (t->nucs[col].allele_cnt != 2)
        return 0;

    this->_logfh << this->_site_index << "\t" << loc->id << "\t" << col << "\t";

    char nuc = '?';

    for (uint p = 0; p < pop_cnt; p++) {
        s = lps->per_pop(p);

        switch(s->nucs[col].p_nuc) {
        case 0:
            nuc = 'N';
            break;
        case 'A':
            switch(s->nucs[col].q_nuc) {
            case 'C':
                nuc = 'M';
                break;
            case 'G':
                nuc = 'R';
                break;
            case 'T':
                nuc = 'W';
                break;
            case 0:
                nuc = 'A';
                break;
            }
            break;
        case 'C':
            switch(s->nucs[col].q_nuc) {
            case 'A':
                nuc = 'M';
                break;
            case 'G':
                nuc = 'S';
                break;
            case 'T':
                nuc = 'Y';
                break;
            case 0:
                nuc = 'C';
                break;
            }
            break;
        case 'G':
            switch(s->nucs[col].q_nuc) {
            case 'A':
                nuc = 'R';
                break;
            case 'C':
                nuc = 'S';
                break;
            case 'T':
                nuc = 'K';
                break;
            case 0:
                nuc = 'G';
                break;
            }
            break;
        case 'T':
            switch(s->nucs[col].q_nuc) {
            case 'A':
                nuc = 'W';
                break;
            case 'C':
                nuc = 'Y';
                break;
            case 'G':
                nuc = 'K';
                break;
            case 0:
                nuc = 'T';
                break;
            }
            break;
        }

        this->_tmpfh << nuc;

        if (p < pop_cnt - 1)
            this->_tmpfh << "\t";

        this->_logfh << this->_mpopi->pops()[p].name << ":" << nuc << ",";

    }
    this->_logfh << "\n";
    this->_site_index++;

    this->_tmpfh << "\n";

    return 0;
}


int
PhylipExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    vector<string> transposed_lines;

    Export::transpose(this->_intmpfh, transposed_lines);

    if (transposed_lines.size() == 0)
        return 0;

    assert(transposed_lines.size() == this->_mpopi->pops().size());

    this->_fh << this->_mpopi->pops().size() << "\t" << this->_site_index << "\n";

    for (size_t line_cnt = 0; line_cnt < transposed_lines.size(); line_cnt++) {
        size_t pos = transposed_lines[line_cnt].find_first_of("\t");
        this->_fh << transposed_lines[line_cnt].substr(0, pos + 1);

        for (size_t i = pos + 1; i < transposed_lines[line_cnt].length(); i++)
            if (transposed_lines[line_cnt][i] != '\t')
                this->_fh << transposed_lines[line_cnt][i];
        this->_fh << "\n";
    }

    this->_fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n";

    return 1;
}

void
PhylipExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    this->_logfh.close();

    remove(this->_tmp_path.c_str());

    this->_fh.close();
    return;
}

int
PhylipVarAllExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // We will write those loci to a Phylip file as defined here:
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    this->_path = out_path + out_prefix + ".all.phylip";
    this->_fh.open(this->_path);
    check_open(this->_fh, this->_path);

    //
    // We will also write a file that allows us to specify each RAD locus as a separate partition
    // for use in phylogenetics programs.
    //
    this->_partition_path = out_path + out_prefix + ".all.partitions.phylip";
    this->_parfh.open(this->_partition_path, ofstream::out);
    check_open(this->_parfh, this->_partition_path);
    
    //
    // Open a log file.
    //
    this->_log_path = this->_path + ".log";
    this->_logfh.open(this->_log_path);

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_logfh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n"
                 << "# Locus ID\tLocusCnt\tSequence position\tLength";
    if (loci_ordered) this->_logfh << "\tChr\tBasepair";
    this->_logfh << "\n";

    cout << "All sites (fixed and variable) will be written in Phylip format to '" << this->_path << "'\n";

    return 0;
}

int
PhylipVarAllExport::write_header()
{
    char   s[id_len];
    size_t len;
    
    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        this->_outstrs[p] += pop.name.substr(0, 10) + "\t";
    }

    sprintf(s, "% 11d", 0);
    this->_fh << this->_mpopi->pops().size() << " " << s << "\n";

    return 0;
}

int
PhylipVarAllExport::write_batch(const vector<LocBin*>& loci)
{
    //
    // We want to write all nucleotides per locus per population in Phylip sequential format. Polymorphic
    // positions will be encoded using IUPAC notation.
    //
    // We will write those loci to a Phylip file as defined here:
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    size_t pop_cnt = this->_mpopi->pops().size();
    char  *seq;

    for (uint k = 0; k < loci.size(); k++) {
        const LocBin*    loc = loci[k];
        const CSLocus*  cloc = loc->cloc;

        // bool include;
        // include = true;
        // for (uint i = 0; i < cloc->snps.size(); i++) {
        //     uint col = cloc->snps[i]->col;

        //     if (t->nucs[col].allele_cnt != 2)
        //         include = false;
        // }
        // if (!include)
        //     continue;

        seq = new char[cloc->len + 1];
        strcpy(seq, cloc->con);

        for (uint j = 0; j < pop_cnt; j++) {
            const LocSum *s = loc->s->per_pop(j);
            
            for (uint snp_index = 0; snp_index < cloc->snps.size(); snp_index++) {
                uint col = cloc->snps[snp_index]->col;

                // if (t->nucs[col].allele_cnt != 2)
                //     continue;
                
                seq[col] = iupac_encode(s->nucs[col].p_nuc, s->nucs[col].q_nuc);
            }

            this->_outstrs[j] += string(seq);
        }
        delete [] seq;

        this->_logfh << cloc->id << "\t" << this->_loc_cnt << "\t" << this->_site_index << "\t" << cloc->len;
        if (loci_ordered) this->_logfh << "\t" << cloc->loc.chr() << "\t" << cloc->sort_bp() + 1;
        this->_logfh << "\n";

        this->_parfh << "DNA, p" << this->_loc_cnt << "=" << this->_site_index << "-" << this->_site_index + cloc->len - 1 << "\n";

        this->_seq_len    += cloc->len;
        this->_site_index += cloc->len;
        this->_loc_cnt++;
    }

    size_t line_len = 120;
    size_t i = 0;
    size_t seql = this->_outstrs.begin()->second.length();
    
    while (i < seql) {
        for (map<int, string>::iterator it = this->_outstrs.begin(); it != this->_outstrs.end(); it++) {
            this->_fh << it->second.substr(i, line_len) << "\n";
        }
        this->_fh << "\n";
        i += line_len;

        if (seql - i < line_len)
            break;
    }

    if (i < seql) {
        for (map<int, string>::iterator it = this->_outstrs.begin(); it != this->_outstrs.end(); it++)
            it->second = it->second.substr(i);
    } else {
        for (map<int, string>::iterator it = this->_outstrs.begin(); it != this->_outstrs.end(); it++)
            it->second = "";
    }
    
    return 0;
}

int
PhylipVarAllExport::post_processing()
{
    //
    // Output any remaining sequence in the buffer.
    //
    if (this->_outstrs.begin()->second.length() > 0) {
        for (map<int, string>::iterator it = this->_outstrs.begin(); it != this->_outstrs.end(); it++)
            this->_fh << it->second << "\n";
    }
    
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the footer.
    //
    this->_fh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n";

    //
    // Rewind the file handle and write the sequence length back over top the first positions.
    //
    this->_fh.seekp(0);

    char   s[id_len];
    sprintf(s, "%11lu", this->_seq_len);
    this->_fh << this->_mpopi->pops().size() << " " << s << "\n";
    
    return 1;
}

void
PhylipVarAllExport::close()
{
    this->_fh.close();
    this->_logfh.close();

    return;
}

int
HzarExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Write a Hybrid Zone Analysis using R (HZAR) file as defined here:
    //    http://cran.r-project.org/web/packages/hzar/hzar.pdf
    //
    this->_path = out_path + out_prefix + ".hzar.csv";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + ".hzar.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "Polymorphic sites in HZAR format will be written to '" << this->_path << "'\n";

    return 0;
}

int
HzarExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    const LocSum   *s;
    const LocTally *t = lps->meta_pop();

    this->_tmpfh << loc->id << "_" << col << ".A"
                 << "," << loc->id << "_" << col << ".B"
                 << "," << loc->id << "_" << col << ".N";

    for (size_t p = 0; p < this->_mpopi->pops().size(); p++) {
        const Pop& pop = this->_mpopi->pops()[p];
        s = lps->per_pop(p);

        this->_tmpfh << "\t" << pop.name;
        
        if (s->nucs[col].num_indv == 0 ||
            s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            //
            // This site contains more than two alleles in this population or was filtered
            // due to a minor allele frequency that is too low.
            //
            this->_tmpfh << "\t" << "0\t0\t0";

        } else {
            if (t->nucs[col].p_allele == s->nucs[col].p_nuc)
                this->_tmpfh << "\t" << s->nucs[col].p << "\t" << 1 - s->nucs[col].p << "\t";
            else
                this->_tmpfh << "\t" << 1 - s->nucs[col].p << "\t" << s->nucs[col].p << "\t";

            this->_tmpfh << s->nucs[col].num_indv * 2;
        }
    }
    this->_tmpfh << "\n";

    return 0;
}

int
HzarExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Reopen the temporary output file to transpose it.
    //
    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    string header;
    map<string, string> transposed_lines;
    vector<string> fields;
    string buf;
    size_t len;
    char   line[max_len];

    uint num_pops = this->_mpopi->pops().size();

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        this->_intmpfh.getline(line, max_len);

        //
        // Check if we read a full line.
        //
        while (this->_intmpfh.fail() && !this->_intmpfh.eof()) {
            buf += line;
            this->_intmpfh.clear();
            this->_intmpfh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        if (!this->_intmpfh.good() || buf.length() == 0)
            break;

        //
        // Break the line up by tabs.
        //
        parse_tsv(buf.c_str(), fields);

        //
        // Each line should have the locus ID, plus four columns per population.
        //
        assert(num_pops * 4 + 1 == fields.size());
        
        header += "," + fields[0];
        
        // Have we seen this population before?
        if (transposed_lines.count(fields[1]) == 0)
            transposed_lines[fields[1]] = "";

        size_t i = 1;
        for (size_t j = 1; j <= num_pops; j++) {
            transposed_lines[fields[i]] += "," + fields[i+1] + "," + fields[i+2] + "," + fields[i+3];
            i += 4;
        }
    } while (!this->_intmpfh.eof());


    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);
    //
    // Output the header line.
    //
    this->_fh << "# Stacks v" << VERSION << "; " << " HZAR v0.2-5; " << date << "\n"
              << "Population" << "," << "Distance" << header << "\n";

    for (size_t p = 0; p < this->_mpopi->pops().size(); p++) {
        const Pop& pop = this->_mpopi->pops()[p];
        this->_fh << pop.name << ",0" << transposed_lines[pop.name] << "\n";
    }
    
    return 1;
}

void
HzarExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();

    remove(this->_tmp_path.c_str());

    this->_fh.close();

    return;
}

int
VcfExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".snps.vcf";
    cout << "SNPs and calls will be written in VCF format to '" << this->_path << "'\n";

    this->_mpopi = mpopi;

    VcfHeader header;
    header.add_std_meta();
    header.add_meta(VcfMeta::predefs::info_loc_strand);
    for(size_t s : this->_mpopi->sample_indexes_orig_order())
        header.add_sample(this->_mpopi->samples()[s].name);

    this->_writer = new VcfWriter(this->_path, move(header));

    return 0;
}

int
VcfExport::write_site(const CSLocus* cloc,
                      const LocPopSum* psum,
                      Datum const*const* d,
                      size_t col,
                      size_t index)
{
    const LocTally* t = psum->meta_pop();

    const char ref = t->nucs[col].p_allele;
    const char alt = t->nucs[col].q_allele;
    char freq_alt[32];
    sprintf(freq_alt, "%0.3f", 1 - t->nucs[col].p_freq);

    VcfRecord rec;
    string id = to_string(cloc->id);
    id += ":";
    id += to_string(col + 1);
    if (loci_ordered) {
        rec.append_chrom(string(cloc->loc.chr()));
        rec.append_pos(cloc->sort_bp(col));
        id += ':';
        id += (cloc->loc.strand == strand_plus ? '+' : '-');
    } else {
        rec.append_chrom(to_string(cloc->id));
        rec.append_pos(col);
    }
    rec.append_id(id);
    rec.append_allele(Nt2(cloc->loc.strand == strand_plus ? ref : reverse(ref)));
    rec.append_allele(Nt2(cloc->loc.strand == strand_plus ? alt : reverse(alt)));
    rec.append_qual();
    rec.append_filters("PASS");
    rec.append_info(string("NS=") + to_string(t->nucs[col].num_indv));
    rec.append_info(string("AF=") + freq_alt);
    rec.append_format("GT");
    rec.append_format("DP");
    rec.append_format("AD");
    rec.append_format("GQ");
    rec.append_format("GL");

    const vector<Nt2> alleles {Nt2(ref), Nt2(alt)};
    for (size_t s : _mpopi->sample_indexes_orig_order()) {
        stringstream sample;

        if (d[s] == NULL || col >= uint(d[s]->len) || d[s]->model[col] == 'U') {
            // Data does not exist.
            sample << "./.";
        } else {
            if (d[s]->model[col] == 'O') {
                assert(d[s]->obshap.size() == 1
                       || (d[s]->obshap.size() == 2 && d[s]->obshap[0][index] == d[s]->obshap[1][index]));
                sample << (d[s]->obshap[0][index] == ref ? "0/0" : "1/1");
            } else {
                assert(d[s]->model[col] == 'E');
                assert(d[s]->obshap.size() == 2 && d[s]->obshap[0][index] != d[s]->obshap[1][index]);
                sample << "0/1";
            }

            // DP.
            sample << ":" << d[s]->snpdata[index].tot_depth;
            // AD.
            if (d[s]->snpdata[index].nt_depths.sum() > 0)
                sample << ":" << d[s]->snpdata[index].nt_depths[Nt2(ref)]
                       << "," << d[s]->snpdata[index].nt_depths[Nt2(alt)];
            else
                sample << ":.";
            // GQ.
            assert(d[s]->snpdata[index].gq != UINT8_MAX);
            sample << ':' << int(d[s]->snpdata[index].gq);
            // GL.
            sample << ':' << VcfRecord::util::fmt_gt_gl(alleles, d[s]->snpdata[index].gtliks);
        }
        rec.append_sample(sample.str());
    }
    this->_writer->write_record(rec);

    return 0;
}

int
VcfHapsExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".haps.vcf";
    cout << "Haplotypes will be written in VCF format to '" << this->_path << "'\n";

    this->_mpopi = mpopi;

    VcfHeader header;
    header.add_std_meta();
    header.add_meta(VcfMeta::predefs::info_loc_strand);
    for(size_t s : this->_mpopi->sample_indexes_orig_order())
        header.add_sample(this->_mpopi->samples()[s].name);

    this->_writer = new VcfWriter(this->_path, move(header));

    return 0;
}

int VcfHapsExport::write_batch(const vector<LocBin*>& loci){
    VcfRecord rec;
    for (const LocBin* locbin : loci) {
        const CSLocus* cloc = locbin->cloc;
        Datum const*const* d = locbin->d;

        if (cloc->snps.empty())
            // Monomorphic locus.
            continue;

        // Tally and sort haplotypes.
        vector<pair<const char*, size_t>> sorted_haps;
        map<const char*, size_t, LessCStrs> hap_indexes;
        tally_complete_haplotypes(d, this->_mpopi->samples().size(), cloc->loc.strand,
                                  sorted_haps, hap_indexes);
        if (sorted_haps.size() < 2)
            continue;

        // Create the record.
        rec.clear();
        if (loci_ordered) {
            rec.append_chrom(string(cloc->loc.chr()));
            rec.append_pos(cloc->loc.bp);
            rec.append_id(to_string(cloc->id) + ":1:" + (cloc->loc.strand == strand_plus ? '+' : '-'));
        } else {
            rec.append_chrom(to_string(cloc->id));
            rec.append_pos(0);
            rec.append_id(".");
        }
        for (size_t i=0; i<sorted_haps.size(); i++) {
            if (cloc->loc.strand == strand_plus) {
                rec.append_allele(string(sorted_haps[i].first));
            } else {
                rec.append_allele(DNASeq4(sorted_haps[i].first).rev_compl().str());
            }
        }
        rec.append_qual();
        rec.append_filters("PASS");
        stringstream info;
        info << "snp_columns=";
        vector<size_t> cols;
        for (const SNP* snp : cloc->snps)
            cols.push_back(snp->col + 1);
        join(cols, ',', info);
        rec.append_info(info.str());
        rec.append_format("GT");

        for(size_t sample : this->_mpopi->sample_indexes_orig_order()) {
            if (d[sample] == NULL) {
                rec.append_sample("./.");
                continue;
            }
            stringstream sample_vcf;
            if (strchr(d[sample]->obshap[0], 'N') != NULL) {
                sample_vcf << "./.";
            } else {
                size_t i0 = hap_indexes.at(d[sample]->obshap[0]);
                size_t i1 = hap_indexes.at(d[sample]->obshap[1]);
                sample_vcf << min(i0, i1) << '/' << max(i0, i1);
            }
            rec.append_sample(sample_vcf.str());
        }
        this->_writer->write_record(rec);
    }

    return 0;
}

int
PlinkExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Write a PLINK file as defined here: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    //

    //
    // Write a markers file containing each marker, the chromosome it falls on,
    // an empty centiMorgan field, and finally its genomic position in basepairs.
    //
    this->_markers_path = out_path + out_prefix + ".plink.map";
    this->_markers_fh.open(this->_markers_path);
    check_open(this->_markers_fh, this->_markers_path);

    //
    // Write the genotypes in a separate file.
    //
    this->_path = out_path + out_prefix + ".plink.ped";
    //
    // Open a temporary file.
    //
    this->_tmpfh.open(this->tmp_path());
    check_open(this->_tmpfh, this->tmp_path());

    cout << "Polymorphic sites in PLINK format will be written to '" << this->_path << "'\n";

    return 0;
}

int
PlinkExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);
    this->_markers_fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    // Population names.
    for (size_t s : this->_mpopi->sample_indexes_orig_order()) {
        const Sample& sample = this->_mpopi->samples()[s];
        this->_tmpfh << '\t' << this->_mpopi->pops()[sample.pop].name;
    }
    this->_tmpfh << '\n';
    // Sample names.
    for (size_t s : this->_mpopi->sample_indexes_orig_order()) {
        const Sample& sample = this->_mpopi->samples()[s];
        this->_tmpfh << '\t' << sample.name;
    }
    this->_tmpfh << '\n';
    // Paternal IDs, maternal IDs, sexes, phenotypes.
    for (size_t i=0; i<4; ++i) {
        for (size_t j=0; j<this->_mpopi->samples().size(); ++j)
            this->_tmpfh << '\t' << 0;
        this->_tmpfh << '\n';
    }
    return 0;
}

int
PlinkExport::write_batch(const vector<LocBin*> &loci)
{
    if (ordered_export) {
        //
        // We need to order the SNPs to take into account overlapping loci.
        //
        vector<const NucTally *> sites;
        OLocTally<NucTally> ord;
        ord.order(sites, loci);

        map<size_t, size_t> key;
        map<size_t, size_t>::iterator key_it = key.begin();

        for (size_t k = 0; k < loci.size(); k++)
            key_it = key.insert(key_it, pair<size_t, size_t>(loci[k]->cloc->id, k));

        for (uint pos = 0; pos < sites.size(); pos++) {
            int loc_id = sites[pos]->loc_id;

            const LocBin* loc = loci[key[loc_id]];
            size_t        col = sites[pos]->col;

            this->_markers_fh << loc->cloc->loc.chr() << "\t"
                              << loc->cloc->id << "_" << col << "\t"
                              << "0\t"
                              << loc->cloc->sort_bp(col) + 1 << "\n";
        }

    } else {
        for (uint k = 0; k < loci.size(); k++) {
            const LocBin*    loc = loci[k];
            const CSLocus*  cloc = loc->cloc;
            const LocTally*    t = loc->s->meta_pop();

            for (uint snp_index = 0; snp_index < cloc->snps.size(); snp_index++) {
                uint col = cloc->snps[snp_index]->col;

                if (t->nucs[col].allele_cnt != 2)
                    continue;

                this->_markers_fh << cloc->loc.chr() << "\t"
                                  << cloc->id << "_" << col << "\t"
                                  << "0\t"
                                  << cloc->sort_bp(col) + 1 << "\n";
            }
        }
    }//

    return OrderableExport::write_batch(loci);
}

int
PlinkExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    const LocSum *s;
    char p_allele, q_allele;

    this->_tmpfh << loc->id << "_" << col;

    for (size_t j : this->_mpopi->sample_indexes_orig_order()) {
        s = lps->per_pop(this->_mpopi->samples()[j].pop);
        if (s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            //
            // This site contains more than two alleles in this population or was filtered
            // due to a minor allele frequency that is too low.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
            //
            // Data does not exist.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j]->model[col] == 'U') {
            //
            // Data exists, but the model call was uncertain.
            //
            this->_tmpfh << "\t" << "0";
        } else {
            //
            // Tally up the nucleotide calls.
            //
            tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

            if (p_allele == 0 && q_allele == 0)
                this->_tmpfh << "\t" << "0";
            else if (p_allele == 0)
                this->_tmpfh << "\t" << q_allele;
            else if (q_allele == 0)
                this->_tmpfh << "\t" << p_allele;
            else
                this->_tmpfh << "\t" << p_allele;
        }
    }

    this->_tmpfh << "\n";

    this->_tmpfh << loc->id << "_" << col;

    for (size_t j : this->_mpopi->sample_indexes_orig_order()) {
        s = lps->per_pop(this->_mpopi->samples()[j].pop);
        if (s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            //
            // This site contains more than two alleles in this population or was filtered
            // due to a minor allele frequency that is too low.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
            //
            // Data does not exist.
            //
            this->_tmpfh << "\t" << "0";
        } else if (d[j]->model[col] == 'U') {
            //
            // Data exists, but the model call was uncertain.
            //
            this->_tmpfh << "\t" << "0";
        } else {
            //
            // Tally up the nucleotide calls.
            //
            tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

            if (p_allele == 0 && q_allele == 0)
                this->_tmpfh << "\t" << "0";
            else if (p_allele == 0)
                this->_tmpfh << "\t" << q_allele;
            else if (q_allele == 0)
                this->_tmpfh << "\t" << p_allele;
            else
                this->_tmpfh << "\t" << q_allele;
        }
    }
    this->_tmpfh << "\n";

    return 0;
}

int
PlinkExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh.open(this->_path.c_str());
    check_open(this->_fh, this->_path);

    //
    // Output the header line.
    //
    this->_fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    ifstream intmpfh (this->tmp_path());
    check_open(intmpfh, this->tmp_path());

    vector<string> transposed_lines;
    Export::transpose(intmpfh, transposed_lines);
    assert(transposed_lines.size() == this->_mpopi->samples().size() + 1);

    //
    // The first line has a list of locus IDs, they are unnecessary in this
    // format, so ignore them.
    //
    auto line = transposed_lines.begin();
    assert(!transposed_lines.empty());
    ++line;
    for (; line != transposed_lines.end(); ++line)
        this->_fh << *line << '\n';

    return true;
}

void
PlinkExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    remove(this->tmp_path().c_str());

    this->_fh.close();
    this->_markers_fh.close();

    return;
}

int
TreemixExport::open(const MetaPopInfo *mpopi)
{
    //
    // Write a TreeMix file (Pickrell and Pritchard, 2012 PLoS Genetics)
    //    https://bitbucket.org/nygcresearch/treemix/wiki/Home
    //
    _mpopi = mpopi;
    _path = out_path + out_prefix + ".treemix";
    cout << "Per-population SNP allele counts will be written to '" << _path << "'\n";
    _writer.open(_path);
    check_open(_writer, _path);

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    _writer << "# Stacks v" << VERSION << "; " << " TreeMix v1.1; " << date << "\n";

    //
    // Output a space-separated list of the populations on the first line.
    //
    bool first = true;
    for (const Pop& p : _mpopi->pops()) {
        if (first)
            first = false;
        else
            _writer << ' ';
        _writer << p.name;

    }
    _writer << '\n';
    return 0;
}

int
TreemixExport::write_site(const CSLocus* cloc,
                      const LocPopSum* psum,
                      Datum const*const* d,
                      size_t col,
                      size_t index)
{
    const LocTally* t = psum->meta_pop();

    bool first = true;
    for (size_t pop=0; pop<_mpopi->pops().size(); ++pop) {
        const LocSum* s = psum->per_pop(pop);
        if (first)
            first = false;
        else
            _writer << ' ';
        if (s->nucs[col].num_indv == 0 ||
            s->nucs[col].incompatible_site ||
            s->nucs[col].filtered_site) {
            _writer << "0,0";
            continue;
        }
        double p_freq = (t->nucs[col].p_allele == s->nucs[col].p_nuc) ?
            s->nucs[col].p :
            1 - s->nucs[col].p;
        size_t allele_cnt = 2 * s->nucs[col].num_indv;
        size_t p_cnt = round(p_freq * allele_cnt);
        size_t q_cnt = allele_cnt - p_cnt;
        _writer << p_cnt << ',' << q_cnt;
    }
    _writer << '\n';
    return 0;
}

int
JoinMapExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Find the 'progeny' population.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "progeny")
            break;
    this->_progeny_index = index;

    this->_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".joinmap.loc";
    
    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".joinmap.loc.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "JoinMap mapping genotypes will be written to '" << this->_path << "'\n";

    //
    // Load the specific marker and genotype dictionaries.
    //
    switch(this->_gproc->type()) {
    case CrossT::cp:
        load_joinmap_cp_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::f2:
        load_mm_f2_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::dh:
        load_mm_dh_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::bc1:
        load_mm_bc_dictionary(this->_marker_map, this->_genotype_map);
        break;
    default:
        cerr << "Unknown mapping cross encountered for output format 'JoinMap'.\n";
        exit(1);
    }

    return 0;
}

int
JoinMapExport::write_batch(const vector<LocBin *> &loci)
{
    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;

    for (uint k = 0; k < loci.size(); k++) {
        const LocBin  *loc  = loci[k];
        const CSLocus *cloc = loc->cloc;

        if (cloc->marker.length() == 0)
            continue;
        
        this->_tmpfh << cloc->id << "\t";

        string joinmap_marker = this->_marker_map[cloc->marker];
        
        if (joinmap_marker == "lmx--")
            this->_tmpfh << "<lmxll>";
        else if (joinmap_marker == "--xnp")
            this->_tmpfh << "<nnxnp>";
        else
            this->_tmpfh << "<" << joinmap_marker << ">";
        
        Datum **d = loc->d;

        for (uint j = fsample; j <= lsample; j++) {
            this->_tmpfh << "\t";

            if (d[j] == NULL) {
                this->_gproc->type() == CrossT::cp ? this->_tmpfh << "--" : this->_tmpfh << "-";

            } else {
                //
                // Assign the program-specific genotype.
                //
                string g = this->_genotype_map[joinmap_marker].count(d[j]->genotype()) ?
                    this->_genotype_map[joinmap_marker][d[j]->genotype()] :
                    this->_genotype_map[joinmap_marker]["--"];
                this->_tmpfh << g;
            }
        }

        this->_tmpfh << "\n";
    }

    return 0;
}

int
JoinMapExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32], date2[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);
    strftime(date2, 32, "%Y%m%d", timeinfo);

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Output the header of the file.
    //
    this->_fh << "# Stacks v" << VERSION << "; JoinMap; " << date << "\n";
    //
    // Output which sample is the male and female in the F0s.
    //
    const vector<Pop>    &pops    = this->_mpopi->pops();
    const vector<Sample> &samples = this->_mpopi->samples();

    size_t index        = this->_mpopi->get_pop_index("parent");
    size_t parent_count = pops[index].n_samples();

    if (parent_count == 1) {
        size_t parent_index = pops[index].first_sample;

        this->_fh << "# Parent: " << samples[parent_index].name << "\n";
    } else {
        size_t index_1 = pops[index].first_sample;
        size_t index_2 = pops[index].last_sample;

        // MetaPopInfo will reorder samples alphabetically, get the original order.
        const vector<size_t> orig_order = this->_mpopi->sample_indexes_orig_order();

        size_t female_index, male_index, i, j;
        for (i = 0; i < orig_order.size(); i++)
            if (orig_order[i] == index_1)
                break;
        for (j = 0; i < orig_order.size(); j++)
            if (orig_order[j] == index_2)
                break;
        if (i < j) {
            female_index = index_1;
            male_index   = index_2;
        } else {
            female_index = index_2;
            male_index   = index_1;
        }
        this->_fh << "# Female [a,b,f,l,m]: " << samples[female_index].name << "\n"
                  << "# Male   [c,d,g,n,p]: " << samples[male_index].name   << "\n";
    }
    this->_fh << "name = populations." << date2 << "\n"
              << "popt = " << this->_gproc->type_str() << "\n"
              << "nloc = " << this->_gproc->mappable_loci() << "\n"
              << "nind = " << this->_mpopi->pops()[this->_progeny_index].n_samples() << "\n\n";

    string buf;
    size_t len;
    char   line[max_len];

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        this->_intmpfh.getline(line, max_len);

        //
        // Check if we read a full line.
        //
        while (this->_intmpfh.fail() && !this->_intmpfh.eof()) {
            buf += line;
            this->_intmpfh.clear();
            this->_intmpfh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        if (!this->_intmpfh.good() || buf.length() == 0)
            break;

        this->_fh << buf << "\n";

    } while (!this->_intmpfh.eof());

    //
    // Write the footer containing the sample names.
    //
    this->_fh << "\nindividual names:\n";
    for (uint j = this->_mpopi->pops()[this->_progeny_index].first_sample; j <= this->_mpopi->pops()[this->_progeny_index].last_sample; j++)
        this->_fh << samples[j].name << "\n";

    return 1;
}

void
JoinMapExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    remove(this->tmp_path().c_str());

    this->_fh.close();

    double progeny_cnt = (double) this->_mpopi->pops()[this->_progeny_index].n_samples();
    
    cerr << "\nJoinMap marker export: "
         << std::setprecision(fieldw)
         << this->_gproc->mappable_loci() << " of " << this->_gproc->n_loci() << " loci were mappable "
         << "(" << (double) this->_gproc->mappable_loci() / (double) this->_gproc->n_loci() * 100 << "%) "
         << "for map type '" << this->_gproc->type_str() << "'; with "
         << this->_gproc->mean_progeny() << " mean mappable progeny per locus "
         << "(" << this->_gproc->mean_progeny() / progeny_cnt * 100 << "%).\n";

    if (this->_gproc->type() == CrossT::cp)
        cerr << this->_gproc->corrected() << " markers were corrected after discovering a missing parental allele.\n";
    
    return;
}

int
OneMapExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Find the 'progeny' population.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "progeny")
            break;
    this->_progeny_index = index;

    this->_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".onemap.txt";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".onemap.txt.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "OneMap mapping genotypes will be written to '" << this->_path << "'\n";

    //
    // Load the specific marker and genotype dictionaries.
    //
    switch(this->_gproc->type()) {
    case CrossT::cp:
        load_onemap_cp_dictionary(this->_marker_map, this->_genotype_map);
        break;
    default:
        cerr << "Unknown mapping cross encountered for output format 'onemap'.\n";
        exit(1);
    }

    return 0;
}

int
OneMapExport::write_batch(const vector<LocBin *> &loci)
{
    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;

    for (uint k = 0; k < loci.size(); k++) {
        const LocBin  *loc  = loci[k];
        const CSLocus *cloc = loc->cloc;

        if (cloc->marker.length() == 0)
            continue;

        string onemap_marker = this->_marker_map[cloc->marker];

        this->_tmpfh << "*" << cloc->id << " " << this->_onemap_marker_types[onemap_marker] << "\t";
        
        Datum **d = loc->d;

        for (uint j = fsample; j <= lsample; j++) {

            if (d[j] == NULL) {
                this->_tmpfh << "-";

            } else {
                //
                // Assign the program-specific genotype.
                //
                string g = this->_genotype_map[onemap_marker].count(d[j]->genotype()) ?
                    this->_genotype_map[onemap_marker][d[j]->genotype()] :
                    this->_genotype_map[onemap_marker]["--"];
                this->_tmpfh << g;

                if (j < lsample)
                    this->_tmpfh << ",";
            }
        }

        this->_tmpfh << "\n";
    }

    return 0;
}

int
OneMapExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Output the header of the file.
    //
    this->_fh << "# Stacks v" << VERSION << "; OneMap; " << date << "\n";
    this->_fh << this->_mpopi->pops()[this->_progeny_index].n_samples() << "\t" << this->_gproc->mappable_loci() << "\n";

    string buf;
    size_t len;
    char   line[max_len];

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        this->_intmpfh.getline(line, max_len);

        //
        // Check if we read a full line.
        //
        while (this->_intmpfh.fail() && !this->_intmpfh.eof()) {
            buf += line;
            this->_intmpfh.clear();
            this->_intmpfh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        if (!this->_intmpfh.good() || buf.length() == 0)
            break;

        this->_fh << buf << "\n";

    } while (!this->_intmpfh.eof());

    return 1;
}

void
OneMapExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    remove(this->tmp_path().c_str());

    this->_fh.close();

    double progeny_cnt = (double) this->_mpopi->pops()[this->_progeny_index].n_samples();
    
    cerr << "\nOneMap marker export: "
         << std::setprecision(fieldw)
         << this->_gproc->mappable_loci() << " of " << this->_gproc->n_loci() << " loci were mappable "
         << "(" << (double) this->_gproc->mappable_loci() / (double) this->_gproc->n_loci() * 100 << "%) "
         << "for map type '" << this->_gproc->type_str() << "'; with "
         << this->_gproc->mean_progeny() << " mean mappable progeny per locus "
         << "(" << this->_gproc->mean_progeny() / progeny_cnt * 100 << "%).\n";

    if (this->_gproc->type() == CrossT::cp)
        cerr << this->_gproc->corrected() << " markers were corrected after discovering a missing parental allele.\n";
    
    return;
}

int
rQTLExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Find the 'progeny' population.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "progeny")
            break;
    this->_progeny_index = index;

    this->_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".rqtl.tsv";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".rqtl.tsv.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "r/QTL mapping genotypes will be written to '" << this->_path << "'\n";

    //
    // Load the specific marker and genotype dictionaries.
    //
    switch(this->_gproc->type()) {
    case CrossT::cp:
        load_joinmap_cp_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::f2:
        load_mm_f2_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::dh:
        load_mm_dh_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::bc1:
        load_mm_bc_dictionary(this->_marker_map, this->_genotype_map);
        break;
    default:
        cerr << "Unknown mapping cross encountered for output format 'rqtl'.\n";
        exit(1);
    }

    return 0;
}

int
rQTLExport::write_header()
{
    this->_tmpfh << "\t\t";

    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;

    for (uint j = fsample; j <= lsample; j++) {
        const Sample& sample = this->_mpopi->samples()[j];
        this->_tmpfh << '\t' << sample.name;
    }
    this->_tmpfh << '\n';

    return 0;
}

int
rQTLExport::write_batch(const vector<LocBin *> &loci)
{
    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;
    
    for (uint k = 0; k < loci.size(); k++) {
        const LocBin  *loc  = loci[k];
        const CSLocus *cloc = loc->cloc;

        if (cloc->marker.length() == 0)
            continue;

        string rqtl_marker = this->_marker_map[cloc->marker];

        this->_tmpfh << cloc->id << "\t" << cloc->loc.chr() << "\t" << cloc->sort_bp();
        
        Datum **d = loc->d;

        for (uint j = fsample; j <= lsample; j++) {
            this->_tmpfh << "\t";
            
            if (d[j] == NULL) {
                this->_gproc->type() == CrossT::cp ? this->_tmpfh << "--" : this->_tmpfh << "-";

            } else {
                //
                // Assign the program-specific genotype.
                //
                string g = this->_genotype_map[rqtl_marker].count(d[j]->genotype()) ?
                    this->_genotype_map[rqtl_marker][d[j]->genotype()] :
                    this->_genotype_map[rqtl_marker]["--"];
                this->_tmpfh << g;
            }
        }

        this->_tmpfh << "\n";
    }

    return 0;
}

int
rQTLExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32], date2[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);
    strftime(date2, 32, "%Y%m%d", timeinfo);

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Output the header of the file.
    //
    this->_fh << "# Stacks v" << VERSION << "; R/QTL; " << date << "\n"
              << "# Exported: populations." << date2 << "\n"
              << "# Map Type: " << this->_gproc->type_str() << "\n"
              << "# Num Loci: " << this->_gproc->mappable_loci() << "\n"
              << "# Num Samples: " << this->_mpopi->pops()[this->_progeny_index].n_samples() << "\n";

    vector<string> transposed_lines;

    Export::transpose(this->_intmpfh, transposed_lines);

    assert(transposed_lines.size() == (this->_mpopi->pops()[this->_progeny_index].n_samples() + 3));

    for (size_t line_cnt = 0; line_cnt < transposed_lines.size(); line_cnt++)
        this->_fh << transposed_lines[line_cnt] << "\n";

    return 1;
}

void
rQTLExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    remove(this->tmp_path().c_str());

    this->_fh.close();

    double progeny_cnt = (double) this->_mpopi->pops()[this->_progeny_index].n_samples();
    
    cerr << "\nR/QTL marker export: "
         << std::setprecision(fieldw)
         << this->_gproc->mappable_loci() << " of " << this->_gproc->n_loci() << " loci were mappable "
         << "(" << (double) this->_gproc->mappable_loci() / (double) this->_gproc->n_loci() * 100 << "%) "
         << "for map type '" << this->_gproc->type_str() << "'; with "
         << this->_gproc->mean_progeny() << " mean mappable progeny per locus "
         << "(" << this->_gproc->mean_progeny() / progeny_cnt * 100 << "%).\n";
    
    return;
}

int
MapMakerExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Find the 'progeny' population.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "progeny")
            break;
    this->_progeny_index = index;

    this->_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".mapmaker.txt";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + "." + this->_gproc->type_str() + ".mapmaker.txt.part";
    this->_tmpfh.open(this->_tmp_path);

    cout << "MapMaker mapping genotypes will be written to '" << this->_path << "'\n";

    //
    // Load the specific marker and genotype dictionaries.
    //
    switch(this->_gproc->type()) {
    case CrossT::f2:
        load_mm_f2_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::dh:
        load_mm_dh_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::bc1:
        load_mm_bc_dictionary(this->_marker_map, this->_genotype_map);
        break;
    default:
        cerr << "Unknown mapping cross encountered for output format 'Mapmaker'.\n";
        exit(1);
    }

    return 0;
}

int
MapMakerExport::write_batch(const vector<LocBin *> &loci)
{
    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;

    for (uint k = 0; k < loci.size(); k++) {
        const LocBin  *loc  = loci[k];
        const CSLocus *cloc = loc->cloc;

        if (cloc->marker.length() == 0)
            continue;

        string mapmaker_marker = this->_marker_map[cloc->marker];
        
        this->_tmpfh << "*" << cloc->id;
        
        Datum **d = loc->d;

        for (uint j = fsample; j <= lsample; j++) {
            this->_tmpfh << " ";

            if (d[j] == NULL) {
                this->_tmpfh << "-";

            } else {
                //
                // Assign the program-specific genotype.
                //
                string g = this->_genotype_map[mapmaker_marker].count(d[j]->genotype()) ?
                    this->_genotype_map[mapmaker_marker][d[j]->genotype()] :
                    this->_genotype_map[mapmaker_marker]["--"];
                this->_tmpfh << g;
            }
        }

        this->_tmpfh << "\n";
    }

    return 0;
}

int
MapMakerExport::post_processing()
{
    //
    // Close the temporary output file.
    //
    this->_tmpfh.close();

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    check_open(this->_intmpfh, this->_tmp_path);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    check_open(this->_fh, this->_path);

    //
    // Output the header of the file.
    //
    this->_fh << "# Stacks v" << VERSION << "; MapMaker; " << date << "\n";

    //
    // Output map type.
    //
    if (this->_gproc->type() == CrossT::f2 )
        this->_fh << "data type f2 intercross\n";
    else if (this->_gproc->type() == CrossT::bc1)
        this->_fh << "data type f2 backcross\n";

    //
    // Output the header: number of individuals, number of markers, number of
    // quantitative traits (none).
    //
    this->_fh << this->_mpopi->pops()[this->_progeny_index].n_samples() << "\t" << this->_gproc->mappable_loci() << "0\n\n";

    string buf;
    size_t len;
    char   line[max_len];

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        this->_intmpfh.getline(line, max_len);

        //
        // Check if we read a full line.
        //
        while (this->_intmpfh.fail() && !this->_intmpfh.eof()) {
            buf += line;
            this->_intmpfh.clear();
            this->_intmpfh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        if (!this->_intmpfh.good() || buf.length() == 0)
            break;

        this->_fh << buf << "\n";

    } while (!this->_intmpfh.eof());

    return 1;
}

void
MapMakerExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();
    remove(this->tmp_path().c_str());

    this->_fh.close();

    double progeny_cnt = (double) this->_mpopi->pops()[this->_progeny_index].n_samples();
    
    cerr << "\nMapMaker marker export: "
         << std::setprecision(fieldw)
         << this->_gproc->mappable_loci() << " of " << this->_gproc->n_loci() << " loci were mappable "
         << "(" << (double) this->_gproc->mappable_loci() / (double) this->_gproc->n_loci() * 100 << "%) "
         << "for map type '" << this->_gproc->type_str() << "'; with "
         << this->_gproc->mean_progeny() << " mean mappable progeny per locus "
         << "(" << this->_gproc->mean_progeny() / progeny_cnt * 100 << "%).\n";
    
    return;
}

/*
int
write_fastphase(map<int, CSLocus *> &catalog,
                PopMap<CSLocus> *pmap,
                PopSum<CSLocus> *psum)
{
    //
    // Write a fastPHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as independent, bi-allelic SNPs. We will write one file per chromosome.
    //
    cout << "Writing population data to fastPHASE files...";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        string file = out_path + out_prefix + "." + it->first + ".fastphase.inp";

        ofstream fh(file.c_str(), ofstream::out);

        if (fh.fail()) {
            cerr << "Error opening fastPHASE file '" << file << "'\n";
            exit(1);
        }

        //
        // Tally up the number of sites
        //
        int  total_sites = 0;
        uint col;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    total_sites++;
            }
        }

        //
        // Output the total number of SNP sites and the number of individuals.
        //
        fh << mpopi.samples().size() << "\n"
           << total_sites    << "\n";

        //
        // We need to determine an ordering that can take into account overlapping RAD sites.
        //
        vector<GenPos> ordered_snps;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    ordered_snps.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
            }
        }
        sort(ordered_snps.begin(), ordered_snps.end());

        //
        // Output the position of each site according to its basepair.
        //
        fh << "P";
        for (uint pos = 0; pos < ordered_snps.size(); pos++) {
            loc = catalog[ordered_snps[pos].id];
            col = loc->snps[ordered_snps[pos].snp_index]->col;
            fh << " " << ordered_snps[pos].bp +1;
        }
        fh << "\n";

        //
        // Output a line of 'S' characters, one per site, indicating that these are SNP markers.
        //
        string snp_markers, gtypes_str;
        snp_markers.assign(total_sites, 'S');
        fh << snp_markers << '\n';

        //
        // Now output each sample name followed by a new line, then all of the genotypes for that sample
        // on two lines.
        //

        char         p_allele, q_allele;
        stringstream gtypes;

        for (size_t p=0; p<mpopi.pops().size(); ++p) {
            const Pop& pop = mpopi.pops()[p];

            for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                //
                // Output all the loci for this sample, printing only the p allele
                //
                fh << mpopi.samples()[j].name << "\n";

                gtypes.str("");
                for (uint pos = 0; pos < ordered_snps.size(); pos++) {
                    loc = catalog[ordered_snps[pos].id];
                    col = loc->snps[ordered_snps[pos].snp_index]->col;

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // If this site is fixed in all populations or has too many alleles don't output it.
                    //
                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        //
                        // This site contains more than two alleles in this population or was filtered
                        // due to a minor allele frequency that is too low.
                        //
                        gtypes << "? ";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        //
                        // Data does not exist.
                        //
                        gtypes << "? ";
                    } else if (d[j]->model[col] == 'U') {
                        //
                        // Data exists, but the model call was uncertain.
                        //
                        gtypes << "? ";
                    } else {
                        //
                        // Tally up the nucleotide calls.
                        //
                        tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            gtypes << "? ";
                        else if (p_allele == 0)
                            gtypes << q_allele << " ";
                        else
                            gtypes << p_allele << " ";
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

                //
                // Output all the loci for this sample again, now for the q allele
                //
                gtypes.str("");
                for (uint pos = 0; pos < ordered_snps.size(); pos++) {
                    loc = catalog[ordered_snps[pos].id];
                    col = loc->snps[ordered_snps[pos].snp_index]->col;


                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        gtypes << "? ";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        gtypes << "? ";

                    } else if (d[j]->model[col] == 'U') {
                        gtypes << "? ";

                    } else {
                        tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            gtypes << "? ";
                        else if (q_allele == 0)
                            gtypes << p_allele << " ";
                        else
                            gtypes << q_allele << " ";
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
            }
        }

        fh.close();
    }

    cout << "done.\n";

    return 0;
}

int
write_phase(map<int, CSLocus *> &catalog,
            PopMap<CSLocus> *pmap,
            PopSum<CSLocus> *psum)
{
    //
    // Write a PHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as mixture of multiple allele, linked RAD sites
    // (SNPs within a single RAD locus are already phased), and bi-allelic SNPs. We
    // will write one file per chromosome.
    //
    cout << "Writing population data to PHASE files...";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        string file = out_path + out_prefix + "." + it->first + ".phase.inp";

        ofstream fh(file.c_str(), ofstream::out);

        if (fh.fail()) {
            cerr << "Error opening PHASE file '" << file << "'\n";
            exit(1);
        }

        //
        // We need to determine an ordering for all legitimate loci/SNPs.
        //
        uint           col;
        vector<GenPos> ordered_loci;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            if (loc->snps.size() == 0) continue;

            //
            // Will we output this locus as a haplotype or as a SNP?
            //
            if (loc->snps.size() > 1) {
                //
                // Check that there aren't too many haplotypes (PHASE has a max of 50).
                //
                if (loc->alleles.size() > 40) continue;

                //
                // Iterate over the population to determine that this subset of the population
                // has data at this locus.
                //
                d = pmap->locus(loc->id);
                for (int j = 0; j < pmap->sample_cnt(); j++) {
                    if (d[j] != NULL &&
                        d[j]->obshap.size() > 0 &&
                        d[j]->obshap.size() <= 2) {
                        //
                        // Data exists, and there are the correct number of haplotypes.
                        //
                        ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(), haplotype));
                        break;
                    }
                }
            } else {
                col = loc->snps[0]->col;

                if (t->nucs[col].allele_cnt == 2)
                    ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(col), snp));
            }
        }
        sort(ordered_loci.begin(), ordered_loci.end());

        //
        // Output the total number of SNP sites and the number of individuals.
        //
        fh << mpopi.samples().size()      << "\n"
           << ordered_loci.size() << "\n";

        //
        // Output the position of each site according to its basepair.
        //
        fh << "P";
        for (uint pos = 0; pos < ordered_loci.size(); pos++)
            fh << " " << ordered_loci[pos].bp +1;
        fh << "\n";

        //
        // Output a line of 'S' characters for SNP markers, 'M' characters for multiallelic haplotypes.
        //
        for (uint pos = 0; pos < ordered_loci.size(); pos++) {
            if (pos > 0) fh << " ";
            fh << (ordered_loci[pos].type == snp ? "S" : "M");
        }
        fh << "\n";

        //
        // Now output each sample name followed by a new line, then all of the genotypes for that sample
        // on two lines.
        //

        string       gtypes_str;
        bool         found;
        char         p_allele, q_allele;
        stringstream gtypes;

        for (size_t p=0; p<mpopi.pops().size(); ++p) {
            const Pop& pop = mpopi.pops()[p];

            for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                //
                // Output all the loci for this sample, printing only the p allele
                //
                fh << mpopi.samples()[j].name << "\n";

                gtypes.str("");
                for (uint pos = 0; pos < ordered_loci.size(); pos++) {
                    loc = catalog[ordered_loci[pos].id];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // Will we output this locus as a haplotype or as a SNP?
                    //
                    if (ordered_loci[pos].type == haplotype) {
                        if (d[j] == NULL) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "-1 ";
                        } else {
                            //
                            // Data exists, output the first haplotype. We will assume the haplotypes are
                            // numbered by their position in the loc->strings vector.
                            //
                            if (d[j]->obshap.size() > 2)  {
                                // cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
                                gtypes << "-1 ";
                            } else {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[0] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Error: Unable to find haplotype " << d[j]->obshap[0] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            }
                        }
                    } else {
                        col = loc->snps[ordered_loci[pos].snp_index]->col;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            //
                            // This site contains more than two alleles in this population or was filtered
                            // due to a minor allele frequency that is too low.
                            //
                            gtypes << "? ";

                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "? ";
                        } else if (d[j]->model[col] == 'U') {
                            //
                            // Data exists, but the model call was uncertain.
                            //
                            gtypes << "? ";
                        } else {
                            //
                            // Tally up the nucleotide calls.
                            //
                            tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                gtypes << "? ";
                            else if (p_allele == 0)
                                gtypes << q_allele << " ";
                            else
                                gtypes << p_allele << " ";
                        }
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

                //
                // Output all the loci for this sample again, now for the q allele
                //
                gtypes.str("");
                for (uint pos = 0; pos < ordered_loci.size(); pos++) {
                    loc = catalog[ordered_loci[pos].id];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // Will we output this locus as a haplotype or as a SNP?
                    //
                    if (ordered_loci[pos].type == haplotype) {
                        if (d[j] == NULL) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "-1 ";
                        } else {
                            //
                            // Data exists, output the second haplotype. We will assume the haplotypes are
                            // numbered by their position in the loc->strings vector.
                            //
                            if (d[j]->obshap.size() > 2)  {
                                // cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
                                gtypes << "-1 ";
                            } else if (d[j]->obshap.size() > 1)  {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[1] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Unable to find haplotype " << d[j]->obshap[1] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            } else {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[0] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Unable to find haplotype " << d[j]->obshap[0] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            }
                        }
                    } else {
                        col = loc->snps[ordered_loci[pos].snp_index]->col;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            gtypes << "? ";

                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            gtypes << "? ";

                        } else if (d[j]->model[col] == 'U') {
                            gtypes << "? ";

                        } else {
                            tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                gtypes << "? ";
                            else if (q_allele == 0)
                                gtypes << p_allele << " ";
                            else
                                gtypes << q_allele << " ";
                        }
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
            }
        }

        fh.close();
    }

    cout << "done.\n";

    return 0;
}
*/

int
tally_observed_haplotypes(const vector<char *> &obshap, int snp_index, char &p_allele, char &q_allele)
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

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
        p_allele = 0;
        q_allele = 0;
        return -1;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    // (The P allele is the first one alphabetically, and the Q allele the second
    // one, if any.)
    //
    p_allele = 0;
    q_allele = 0;

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

    return 0;
}
