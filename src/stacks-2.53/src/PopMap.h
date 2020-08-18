// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __POPMAP_H__
#define __POPMAP_H__

#include <exception>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <numeric>
using std::accumulate;
#include <algorithm>
#include <utility>

#include "stacks.h"
#include "locus.h"
#include "aln_utils.h"
#include "MetaPopInfo.h"
#include "Vcf.h"

class Datum {
public:
    struct SNPData {
        size_t tot_depth;
        Counts<Nt2> nt_depths;
        uint8_t gq;
        GtLiks gtliks;
        SNPData() : tot_depth(0), gq(UINT8_MAX) {}
    };

    int             id;            // Stack ID
    int             merge_partner; // Stack ID of merged datum, if this datum was merged/phased from two, overlapping datums.
    int             len;           // Length of locus
    char           *model;         // String representing SNP model output for each nucleotide at this locus.
    char           *gtype;         // Mapping genotype
    char           *cigar;         // CIGAR string describing how the datum aligns to the catalog locus.
    vector<char *>  obshap;        // Observed Haplotypes
    vector<SNPData> snpdata;

    Datum()  {
        this->id            = -1;
        this->gtype         = NULL;
        this->model         = NULL;
        this->cigar         = NULL;
        this->len           = 0;
        this->merge_partner = 0;
    }
    ~Datum() {
        for (uint i = 0; i < this->obshap.size(); i++)
            delete [] this->obshap[i];
        delete [] this->gtype;
        delete [] this->model;
        delete [] this->cigar;
    }
    
    void genotype(string g) {
        if (g.length() == 0) return;
        if (this->gtype != NULL) delete [] this->gtype;
        this->gtype = new char [g.length() + 1];
        strcpy(this->gtype, g.c_str());
    }
    string genotype() { return string(this->gtype); }
};

template<class LocusT=Locus>
class PopMap {
    const MetaPopInfo& metapopinfo;
    int      num_loci;
    set<pair<int, int> > blacklist;
    Datum ***data;
    map<int, int> locus_order;  // LocusID => ArrayIndex; map catalog IDs to their first dimension
                                // position in the Datum array.
    map<int, int> rev_locus_order;
    map<string, vector<LocusT *> > ordered_loci_; // Loci ordered by genomic position

public:
    PopMap(const MetaPopInfo& mpopi, int n_loci);
    ~PopMap();

    // Populates the Locus & PopMap based on Stacks (v2) files.
    static void populate_internal(
        LocusT* cloc,
        Datum** locdata,
        const Seq& fasta_record,
        const vector<VcfRecord>& vcf_records,
        const VcfHeader& vcf_header,
        const MetaPopInfo* mpopi,
        const vector<size_t>& sample_vcf_i_to_mpopi_i
    );
    // Populates the Locus & PopMap based on external VCF files. Returns false
    // and prints warnings if the parsing of the record fails.
    static bool populate_external(
        LocusT* cloc,
        Datum** locdata,
        int cloc_id,
        const VcfRecord& vcf_record,
        const VcfHeader& vcf_header,
        const MetaPopInfo* mpopi,
        const vector<size_t>& sample_vcf_i_to_mpopi_i
    );

    //
    // Obtain indexes, etc.
    //
    int loci_cnt() const { return this->num_loci; }
    int locus_index(int id) const {return locus_order.at(id);}
    int rev_locus_index(int index) const {try {return rev_locus_order.at(index);} catch(exception&) {return -1;}}

    int sample_cnt() const { return metapopinfo.samples().size(); }
    int sample_index(int id) const {try {return metapopinfo.get_sample_index(id);} catch (exception&) {return -1;}}
    int rev_sample_index(int index) const {return metapopinfo.samples().at(index).id;}

    bool blacklisted(int loc_id, int sample_id) const {return blacklist.count({sample_id, loc_id});}

    //
    // Access the Datums.
    //
    Datum **locus(int id) { return this->data[this->locus_order[id]]; }
    Datum  *datum(int loc_id, int sample_id) { return this->locus(loc_id)[metapopinfo.get_sample_index(sample_id)]; }
};

template<class LocusT>
PopMap<LocusT>::PopMap(const MetaPopInfo& mpopi, int num_loci): metapopinfo(mpopi)
{
    this->data = new Datum **[num_loci];

    for (int i = 0; i < num_loci; i++) {
        this->data[i] = new Datum *[metapopinfo.samples().size()];

        for (size_t j = 0; j < metapopinfo.samples().size(); j++)
            this->data[i][j] = NULL;
    }

    this->num_loci    = num_loci;
}

template<class LocusT>
PopMap<LocusT>::~PopMap() {
    for (int i = 0; i < this->num_loci; i++) {
        for (size_t j = 0; j < metapopinfo.samples().size(); j++)
            delete this->data[i][j];
        delete [] this->data[i];
    }
    delete [] this->data;
}

template<class LocusT>
void PopMap<LocusT>::populate_internal(
    LocusT* cloc,
    Datum** locdata,
    const Seq& fasta_record,
    const vector<VcfRecord>& vcf_records,
    const VcfHeader& vcf_header,
    const MetaPopInfo* mpopi,
    const vector<size_t>& sample_vcf_i_to_mpopi_i
) { try {
    assert(fasta_record.id != NULL);
    assert(fasta_record.seq != NULL);
    assert(!fasta_record.comment.empty());
    assert(!vcf_records.empty());
    assert(sample_vcf_i_to_mpopi_i.size() == vcf_header.samples().size());

    // Parse the FASTA record.
    // ==========
    // Locus ID.
    cloc->id = is_integer(fasta_record.id);
    if (cloc->id < 0) {
        cerr << "Error: Unable to parse locus ID.\n";
        throw exception();
    }
    // Consensus sequence.
    cloc->add_consensus(fasta_record.seq);
    if (cloc->len < vcf_records.back().pos() + 1) {
        cerr << "Error: Bad consensus length.\n";
        throw exception();
    }
    // Locus position:
    // If the analysis is reference-based, there will be a ' pos=...' field on
    // the fasta ID line, placed as a comment adjacent to the locus ID.
    // Overlap: if this data is de novo, there will potentially be an overlap
    // length between the single and paired-end contigs.
    assert(cloc->loc.empty());
    const char *p, *q;
    p = fasta_record.comment.c_str();
    do {
        q = p;
        while (*q != '\0' && *q != ' ' && *q != '\t')
            ++q;
        if (strncmp(p, "pos=", 4) == 0) {
            p += 4;
            cloc->loc = PhyLoc(string(p, q));
        } else if (strncmp(p, "contig=overlapped:", 18) == 0) {
            p += 18;
            cloc->overlap = is_integer(p);
            cloc->pe_ctg = true;
        } else if (strncmp(p, "contig=separate", 15) == 0) {
            cloc->pe_ctg = true;
        }
        p = q;
        while(*p != '\0' && (*p == ' ' || *p == '\t'))
            ++p;
    } while (*p != '\0');
    if (cloc->loc.empty())
        cloc->loc = PhyLoc("", 0, strand_plus);

    // Parse the VCF records.
    // ==========

    // CSLocus.
    // ----------
    for (auto& rec : vcf_records) {
        if (!rec.is_monomorphic()) {
            SNP* snp = new SNP();
            cloc->snps.push_back(snp);
            snp->type = snp_type_het;
            snp->col = rec.pos();
            auto a = rec.begin_alleles();
            snp->rank_1 = **a;
            ++a;
            assert(a != rec.end_alleles());
            snp->rank_2 = **a;
            ++a;
            if (a != rec.end_alleles()) {
                snp->rank_3 = **a;
                ++a;
                if (a != rec.end_alleles())
                    snp->rank_4 = **a;
            }
        }
    }

    // Genotypes.
    // ----------
    size_t snp_i = 0;
    --snp_i;
    vector<Nt2> rec_alleles;
    size_t dp_index = SIZE_MAX;
    size_t ad_index = SIZE_MAX;
    size_t gq_index = SIZE_MAX;
    size_t gl_index = SIZE_MAX;
    for (auto& rec : vcf_records)
    { try {
        size_t col = rec.pos();
        bool snp_rec = !rec.is_monomorphic();
        if (snp_rec) {
            if (!rec.is_snp()) {
                cerr << "Error: Not a SNP.\n";
                throw exception();
            }
            ++snp_i;
            rec_alleles.clear();
            for(auto a = rec.begin_alleles(); a != rec.end_alleles(); ++a)
                rec_alleles.push_back(Nt2(**a));
            if (rec.index_of_gt_subfield("PS") != 1)
                throw exception();
            dp_index = rec.index_of_gt_subfield("DP");
            ad_index = rec.index_of_gt_subfield("AD");
            gq_index = rec.index_of_gt_subfield("GQ");
            gl_index = rec.index_of_gt_subfield("GL");
        } else {
            assert(rec.count_formats() == 1 && strcmp(rec.format0(),"DP")==0);
        }
        VcfRecord::iterator gt_itr = rec.begin_samples();
        for (size_t sample_vcf_i = 0;
            sample_vcf_i < vcf_header.samples().size();
            ++gt_itr, ++sample_vcf_i
        ) { try {
            assert(gt_itr != rec.end_samples());
            // Check that the sample is present in the population map.
            size_t sample_mpopi_i = sample_vcf_i_to_mpopi_i[sample_vcf_i];
            if (sample_mpopi_i == SIZE_MAX)
                continue;
            else
                assert(sample_mpopi_i < mpopi->samples().size());
            // Check if the sample has data.
            const char* gt_str = *gt_itr;
            if (gt_str[0] == '.')
                continue;
            Datum* d = locdata[sample_mpopi_i];
            if (snp_rec) {
                // Check that this isn't an unphased HET.
                // N.B.: GStacks writes a single phase set per sample.
                pair<long,long> gt = VcfRecord::util::parse_gt_gt(gt_str);
                assert(gt.first >= 0);
                if (gt.first != gt.second) {
                    const char* ps = strchr(gt_str, ':');
                    if (ps == NULL)
                        throw exception();
                    ++ps;
                    if (ps[0] == '.')
                        continue;
                }
            }
            // This sample has data for this locus for this column.
            if (d == NULL) {
                ++cloc->cnt;
                d = new Datum();
                locdata[sample_mpopi_i] = d;
                d->id = cloc->id;
                d->len = cloc->len;
                d->model = new char[cloc->len+1];
                memset(d->model, 'U', cloc->len);
                d->model[cloc->len] = '\0';
                if (!cloc->snps.empty()) {
                    d->snpdata = vector<Datum::SNPData>(cloc->snps.size());
                    for (size_t i=0; i<2; ++i) {
                        char* hap = new char[cloc->snps.size()+1];
                        d->obshap.push_back(hap);
                        memset(hap, 'N', cloc->snps.size());
                        hap[cloc->snps.size()] = '\0';
                    }
                } else {
                    d->obshap.push_back(new char[10]);
                    strncpy(d->obshap[0], "consensus", 10);
                }
            }
            // Handle the non-SNP-record case.
            if (!snp_rec) {
                d->model[col] = 'O';
                continue;
            }
            // SNP column.
            pair<long,long> gt = VcfRecord::util::parse_gt_gt(gt_str);
            d->model[col] = (gt.first == gt.second ? 'O' : 'E');
            d->obshap[0][snp_i] = char(rec_alleles.at(gt.first));
            d->obshap[1][snp_i] = char(rec_alleles.at(gt.second));
            // Record additional information on the call.
            Datum::SNPData& s = d->snpdata[snp_i];
            s.tot_depth = VcfRecord::util::parse_gt_dp(gt_str, dp_index);
            s.nt_depths = VcfRecord::util::parse_gt_ad(gt_str, ad_index, rec_alleles);
            s.gq = VcfRecord::util::parse_gt_gq(gt_str, gq_index);
            s.gtliks = VcfRecord::util::parse_gt_gl(gt_str, gl_index, rec_alleles);
        } catch (exception&) {
            cerr << "Error: At the " << (sample_vcf_i+1) << "th sample.\n";
            throw;
        }}
    } catch (exception&) {
        cerr << "Error: In record '" << rec.chrom() << ":" << rec.pos()+1 << "'.\n";
        throw;
    }}

    // Finalize the CSLocus (??).
    // ==========
    if (!cloc->snps.empty()) {
        string hap;
        for (size_t i = 0; i < mpopi->samples().size(); ++i) {
            Datum* d = locdata[i];
            if (d == NULL)
                continue;
            assert(d->obshap.size() == 2);
            if (strchr(d->obshap[0], 'N') == NULL) {
                assert(strchr(d->obshap[1], 'N') == NULL);
                hap = d->obshap[0];
                ++cloc->alleles[move(hap)];
                hap = d->obshap[1];
                ++cloc->alleles[move(hap)];
            }
        }
        cloc->populate_alleles();
    }

} catch (exception&) {
    cerr << "Error: Locus " << fasta_record.id << "\n"
         << "Error: Bad GStacks files.\n";
    throw;
}}

template<class LocusT> bool
PopMap<LocusT>::populate_external(
    LocusT* cloc,
    Datum** locdata,
    int cloc_id,
    const VcfRecord& vcf_record,
    const VcfHeader& vcf_header,
    const MetaPopInfo* mpopi,
    const vector<size_t>& sample_vcf_i_to_mpopi_i
) { try {
    assert(vcf_record.is_snp());
    assert(vcf_record.count_alleles() == 2);
    assert(sample_vcf_i_to_mpopi_i.size() == vcf_header.samples().size());
    // We ignore the '*' allele & treat samples that have it as missing.
    long upstream_del_allele = -1;
    vector<Nt2> rec_alleles;
    long allele_i = 0;
try {
    for (auto a = vcf_record.begin_alleles();
        a != vcf_record.end_alleles();
        ++a, ++allele_i
    ) {
        if (**a == '*') {
            if (upstream_del_allele != -1)
                throw exception();
            upstream_del_allele = allele_i;
            rec_alleles.push_back(Nt2::a);
        } else {
            if (!Nt4(**a).is_acgt())
                throw exception();
            rec_alleles.push_back(Nt2(**a));
        }
    }
    if (rec_alleles.size() > 4 + (upstream_del_allele != -1 ? 1 : 0))
        throw exception();
} catch (exception&) {
    cerr << "Warning: Illegal alleles.\n";
    throw;
}
    // CSLocus.
    // ==========
    cloc->id = cloc_id;
    cloc->len = 1;
    cloc->con = new char[2];
    cloc->con[0] = *vcf_record.allele0();
    cloc->con[1] = '\0';
    cloc->loc.set(vcf_record.chrom(), vcf_record.pos(), strand_plus);
    // Locus::snps
    SNP* snp = new SNP();
    cloc->snps.push_back(snp);
    snp->type = snp_type_het;
    snp->col = 0;
    // Locus::rankX
    allele_i = 0;
    if (allele_i == upstream_del_allele)
        ++allele_i;
    snp->rank_1 = char(rec_alleles[allele_i]);
    ++allele_i;
    if (allele_i == upstream_del_allele)
        ++allele_i;
    assert(allele_i != long(rec_alleles.size()));
    snp->rank_2 = char(rec_alleles[allele_i]);
    ++allele_i;
    if (allele_i == upstream_del_allele)
        ++allele_i;
    if (allele_i != long(rec_alleles.size())) {
        snp->rank_3 = char(rec_alleles[allele_i]);
        ++allele_i;
        if (allele_i == upstream_del_allele)
            ++allele_i;
        if (allele_i != long(rec_alleles.size()))
            snp->rank_4 = char(rec_alleles[allele_i]);
    }

    // Genotypes.
    // ==========
    size_t dp_index = vcf_record.index_of_gt_subfield("DP");
    size_t ad_index = vcf_record.index_of_gt_subfield("AD");
    size_t gq_index = vcf_record.index_of_gt_subfield("GQ");
    size_t gl_index = vcf_record.index_of_gt_subfield("GL");
    VcfRecord::iterator gt_itr = vcf_record.begin_samples();
    for (size_t sample_vcf_i = 0;
        sample_vcf_i < vcf_header.samples().size();
        ++gt_itr, ++sample_vcf_i
    ) { try {
        assert(gt_itr != vcf_record.end_samples());
        // Check that the sample is present in the population map.
        size_t sample_mpopi_i = sample_vcf_i_to_mpopi_i[sample_vcf_i];
        if (sample_mpopi_i == SIZE_MAX)
            continue;
        else
            assert(sample_mpopi_i < mpopi->samples().size());
        // Check if the sample has data.
        const char* gt_str = *gt_itr;
        pair<long,long> gt = VcfRecord::util::parse_gt_gt(gt_str);
        if (gt.first == -1)
            continue;
        else if (gt.first == upstream_del_allele || gt.second == upstream_del_allele)
            continue;
        // Fill the datum.
        assert(locdata[sample_mpopi_i] == NULL);
        ++cloc->cnt;
        Datum* d = new Datum();
        locdata[sample_mpopi_i] = d;
        d->id = cloc->id;
        d->len = 1;
        d->model = new char[2];
        d->model[0] = (gt.first == gt.second ? 'O' : 'E');
        d->model[1] = '\0';
        for (size_t i=0; i<2; ++i) {
            d->obshap.push_back(new char[2]);
            d->obshap.back()[0] = 'N';
            d->obshap.back()[1] = '\0';
        }
        d->obshap[0][0] = char(rec_alleles.at(gt.first));
        d->obshap[1][0] = char(rec_alleles.at(gt.second));
        d->snpdata = vector<Datum::SNPData>(1);
        Datum::SNPData& s = d->snpdata[0];
        if (dp_index != SIZE_MAX)
            s.tot_depth = VcfRecord::util::parse_gt_dp(gt_str, dp_index);
        if (ad_index != SIZE_MAX)
            s.nt_depths = VcfRecord::util::parse_gt_ad(gt_str, ad_index, rec_alleles);
        if (gq_index != SIZE_MAX)
            s.gq = VcfRecord::util::parse_gt_gq(gt_str, gq_index);
        if (gl_index != SIZE_MAX)
            s.gtliks = VcfRecord::util::parse_gt_gl(gt_str, gl_index, rec_alleles);
    } catch (exception&) {
        cerr << "Warning: Malformed sample field '" << *gt_itr << "'.\n";
        throw;
    }}

    // Finally, "populate alleles".
    for (size_t i = 0; i < mpopi->samples().size(); ++i) {
        Datum* d = locdata[i];
        if (d == NULL)
            continue;
        assert(d->obshap.size() == 2);
        assert(strchr(d->obshap[0], 'N') == NULL);
        ++cloc->alleles[d->obshap[0]];
        ++cloc->alleles[d->obshap[1]];
    }
    cloc->populate_alleles();
    return true;
} catch (exception&) {
    cerr << "Warning: Discarding VCF SNP record '"
         << vcf_record.chrom() << ":" << (vcf_record.pos()+1) << "'.\n";
    return false;
}}

#endif // __POPMAP_H__
