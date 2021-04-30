// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include "constants.h"
#include "stacks.h"
#include "MetaPopInfo.h"
#include "Alignment.h"
#include "aln_utils.h"

class Match {
 public:
    uint        cat_id;
    allele_type cat_type;
    allele_type query_type;
    string      cigar;
    uint        dist;

    Match() {
        this->cat_id = 0;
        this->dist   = 0;
    };
    Match(const Match &other) {
        this->cat_id     = other.cat_id;
        this->cat_type   = other.cat_type;
        this->query_type = other.query_type;
        this->cigar      = other.cigar;
        this->dist       = other.dist;
    };
};

class Locus;
int  adjust_snps_for_gaps(Cigar&, Locus*);
int  adjust_and_add_snps_for_gaps(Cigar&, Locus*);
int  remove_snps_from_gaps(Cigar&, Locus*);

class Locus {
 public:
    int          id; // Locus ID
    int   sample_id; // Sample ID
    int       depth; // Stack depth
    char       *con; // Consensus sequence
    char     *model; // Model calls for each nucleotide
    uint        len; // Sequence length

    //
    // Flags
    //
    bool blacklisted;
    bool deleveraged;
    bool lumberjackstack;

    vector<char *>          comp;   // Raw components in this stack.
    vector<char *>         reads;   // Sequence reads contributing to this locus.
    vector<uint>        comp_cnt;   // Counter for internal stacks merged into this locus.
    vector<read_type>  comp_type;   // Read types for reads contributing to this locus.
    PhyLoc                   loc;   // Physical genome location of this stack.
    vector<SNP *>           snps;   // Single Nucleotide Polymorphisms in this stack.
    map<string, int>     alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
    vector<pair<allele_type, string>> strings; // Strings for matching (representing the various allele combinations)

    Locus()  {
        id              = 0;
        sample_id       = 0;
        depth           = 0;
        model           = NULL;
        con             = NULL;
        len             = 0;
        blacklisted     = false;
        deleveraged     = false;
        lumberjackstack = false;
    }
    Locus(const Locus &other);
    virtual ~Locus() {
        delete [] con;
        delete [] model;
        for (uint i = 0; i < snps.size(); i++)
            delete snps[i];
        for (uint i = 0; i < comp.size(); i++)
            delete [] comp[i];
        for (uint i = 0; i < reads.size(); i++)
            delete [] reads[i];
    }
    uint sort_bp() const {return this->loc.bp;}
    uint sort_bp(uint col0) const {return this->loc.bp + (this->loc.strand==strand_plus ? col0 : -col0);}
    int snp_index(uint) const;
    int add_consensus(const char *);
    int add_model(const char *);
    virtual int populate_alleles();
};

//
// Query Locus Class
//
class QLocus : public Locus {
 public:
    vector<Match *> matches;   // Matching tags found for the catalog.

    QLocus(): Locus() {}
    QLocus(const QLocus &other);
    ~QLocus();

    int add_match(int, allele_type, allele_type, int, string);
    int add_match(int, allele_type, allele_type, int);
    int add_match(int, allele_type);
    int clear_matches();
};

//
// Catalog Locus Class, for use in cstacks, contains catalog loci and records the
// constiuent tags this locus was built from.
//
class CLocus : public Locus {
 public:
    vector<pair<int, int> > sources;   // Sample/ID pairs for the sources contributing to this catalog entry
    uint match_cnt;

    CLocus() : Locus() {
        this->match_cnt = 0;
    };

    int merge_snps(Locus *);
    int reduce_alleles(set<string> &);
};

//
// Catalog Summary Locus Class; used in genotypes and populations, records a catalog
// locus with summary information derived from individuals in the population.
//
class CSLocus : public Locus {
public:
    CSLocus() : Locus() {
        this->chisq          = 1.0;
        this->cnt            = 0;
        this->gcnt           = 0;
        this->pe_ctg         = false;
        this->overlap        = 0;
    };
    string   annotation;
    string   marker;
    string   uncor_marker;
    map<string, int> hap_cnts; // Counts of each observed haplotype for this locus in the population.
    bool     pe_ctg;           // Was a paired-end contig constructed for this locus?
    uint16_t overlap;          // Size of overlap between single and paired-end contigs for this locus.
    uint16_t cnt;              // Number of samples containing data for this locus.
    uint16_t gcnt;             // Number of progeny containing a valid mapping genotype.
    double   chisq;            // Chi squared p-value testing the null hypothesis of no segregation distortion.
};

// SRead: a Read belonging to a Sample.
struct SRead : Read {
    size_t sample; // index in MetaPopInfo::samples_
    SRead() : Read(), sample(SIZE_MAX) {}
    SRead(Read&& r, size_t spl) : Read(move(r)), sample(spl) {}
};

struct SAlnRead : AlnRead {
    size_t sample; // index in MetaPopInfo::samples_
    SAlnRead() : AlnRead(), sample(SIZE_MAX) {}
    SAlnRead(AlnRead&& r, size_t spl) : AlnRead(move(r)), sample(spl) {}
    SAlnRead(Read&& r, Cigar&& c, size_t spl) : AlnRead(move(r), move(c)), sample(spl) {}
};

class CLocReadSet {
    const MetaPopInfo* mpopi_;
    size_t bam_i_; // BAM target index (SIZE_MAX if unavailable).
    int id_; // Catalog locus ID
    PhyLoc aln_pos_;
    vector<SRead> reads_; // Forward reads. Order is arbitrary.
    vector<SRead> pe_reads_; // Paired-end reads. Order and size are arbitrary.

public:
    CLocReadSet(const MetaPopInfo& mpopi)
        : mpopi_(&mpopi), bam_i_(SIZE_MAX), id_(-1), aln_pos_(), reads_(), pe_reads_()
        {}

    const MetaPopInfo& mpopi() const {return *mpopi_;}
    size_t bam_i() const {return bam_i_;}
    int id() const {return id_;}
    const PhyLoc& pos() const {return aln_pos_;}
    const vector<SRead>& reads() const {return reads_;}
          vector<SRead>& reads()       {return reads_;}
    const vector<SRead>& pe_reads() const {return pe_reads_;}
          vector<SRead>& pe_reads()       {return pe_reads_;}

    void clear() {bam_i_=SIZE_MAX; id_=-1; reads_.clear(); pe_reads_.clear();}
    void bam_i(size_t i) {bam_i_ = i;}
    void id(int id) {id_ = id;}
    void pos(const PhyLoc& p) {aln_pos_ = p;}
    void add(SRead&& r) {reads_.push_back(move(r));}
    void add_pe(SRead&& r) {pe_reads_.push_back(move(r));}
    size_t n_samples() const;
};

class CLocAlnSet {
    int id_; // Catalog locus ID
    PhyLoc aln_pos_;
    DNASeq4 ref_;
    const MetaPopInfo* mpopi_;
    vector<SAlnRead> reads_;
    vector<vector<size_t>> reads_per_sample_; // `at(sample)` is a vector of indexes in `reads_`.

public:
    CLocAlnSet()
        : id_(-1), mpopi_(NULL)
        {}
    CLocAlnSet(CLocAlnSet&&) = default;
    CLocAlnSet& operator= (CLocAlnSet&&) = default;

    const MetaPopInfo& mpopi() const {return *mpopi_;}
    int id() const {return id_;}
    const PhyLoc& pos() const {return aln_pos_;}
    const DNASeq4& ref() const {return ref_;}
    const vector<SAlnRead>& reads() const {return reads_;}
    const vector<size_t>& sample_reads(size_t sample) const {return reads_per_sample_.at(sample);}

    void clear();
    void reinit(int id, const PhyLoc& aln_pos, const MetaPopInfo* mpopi);
    void ref(DNASeq4&& s) {ref_ = move(s);}
    void add(SAlnRead&& r);

    void recompute_consensus();
    void sort_by_read_name();
    void sort_by_alignment_offset();
    void hard_clip_right_Ns();
    void merge_paired_reads();
    void remove_unmerged_reads(ostream* log = NULL);
    void remove_pcr_duplicates(vector<size_t>* clone_size_distrib = NULL, ostream* log = NULL);
    size_t n_samples() const {size_t n=0; for(auto& reads : reads_per_sample_) if (!reads.empty()) ++n; return n;}

    friend ostream& operator<< (ostream& os, const CLocAlnSet& loc);

    static CLocAlnSet juxtapose(CLocAlnSet&& left, CLocAlnSet&& right, long offset);

    //
    // Class to iterate over sites.
    //
    class site_iterator {

        const CLocAlnSet& loc_aln_;
        DNASeq4::iterator ref_it_;
        DNASeq4::iterator ref_past_;
        vector<Alignment::iterator> its_;
        size_t col_;

    public:
        // Iteration methods.
        site_iterator(const CLocAlnSet& loc_aln)
                : loc_aln_(loc_aln),
                ref_it_(loc_aln.ref().begin()),
                ref_past_(loc_aln.ref().end()),
                col_(0)
                {
            its_.reserve(loc_aln.reads().size());
            for (const SAlnRead& r: loc_aln.reads())
                its_.push_back(Alignment::iterator(r.aln()));
        }
        operator bool () const {return ref_it_ != ref_past_;}
        site_iterator& operator++ ();

        // Site interface.
        Nt4 ref_nt() const {return *ref_it_;} // Get the contig nt4.
        void counts(Counts<Nt4>& counts) const; // Get the nt counts across all samples.
        void counts(Counts<Nt4>& counts, size_t sample) const; // Get the nt counts for a given sample.
        SiteCounts counts() const;
        size_t col() const {return col_;}

        const MetaPopInfo& mpopi() const {return loc_aln_.mpopi();}
    };
};

//
// ==================
// Inline definitions
// ==================
//

inline
void CLocAlnSet::add(SAlnRead&& r) {
    assert(cigar_length_ref(r.cigar) == ref_.length());
    reads_per_sample_.at(r.sample).push_back(reads_.size());
    reads_.push_back(move(r));
}

inline
CLocAlnSet::site_iterator& CLocAlnSet::site_iterator::operator++ () {
    ++ref_it_;
    for (auto& it: its_)
        ++it;

    #ifdef DEBUG
    if (! (ref_it_ != ref_past_)) {
        // Make sure we've reached the end of the alignment for every read.
        for (auto& it: its_)
            assert(!bool(it));
    }
    #endif

    ++col_;
    return *this;
}

inline
void CLocAlnSet::site_iterator::counts(Counts<Nt4>& counts) const {
    counts.clear();
    for (auto& read: its_)
        counts.increment(*read);
}

inline
void CLocAlnSet::site_iterator::counts(Counts<Nt4>& counts, size_t sample) const {
    counts.clear();
    for (size_t read_i : loc_aln_.sample_reads(sample))
        counts.increment(*its_[read_i]);
}

inline
SiteCounts CLocAlnSet::site_iterator::counts() const {
    SiteCounts cnts;
    cnts.mpopi = &mpopi();

    cnts.samples.reserve(mpopi().samples().size());
    Counts<Nt4> tmp;
    for (size_t sample=0; sample<mpopi().samples().size(); ++sample) {
        counts(tmp, sample); //this->counts()
        cnts.samples.push_back(Counts<Nt2>(tmp));
        cnts.tot += cnts.samples.back();
    }
    return cnts;
}

#endif // __LOCUS_H__
