// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-OA
//
// Copyright 2013-2018, Julian Catchen <jcatchen@illinois.edu>
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
// locus.cc -- routines for the Locus class and its derivatives.
//
#include "locus.h"

#include "utils.h"

Locus::Locus(const Locus &templ)
{
    this->id              = templ.id;
    this->sample_id       = templ.sample_id;
    this->depth           = templ.depth;
    this->len             = templ.len;
    this->blacklisted     = templ.blacklisted;
    this->deleveraged     = templ.deleveraged;
    this->lumberjackstack = templ.lumberjackstack;

    this->model = NULL;
    this->con   = NULL;

    if (templ.model != NULL) {
        this->model = new char [templ.len + 1];
        strcpy(this->model, templ.model);
    }
    if (templ.con != NULL) {
        this->con = new char [templ.len + 1];
        strcpy(this->con, templ.con);
    }

    for (uint i = 0; i < templ.comp.size(); i++) {
        char *c = new char [strlen(templ.comp[i]) + 1];
        strcpy(c, templ.comp[i]);
        this->comp.push_back(c);
    }
    for (uint i = 0; i < templ.reads.size(); i++) {
        char *c = new char [strlen(templ.reads[i]) + 1];
        strcpy(c, templ.reads[i]);
        this->reads.push_back(c);
    }
    for (uint i = 0; i < templ.snps.size(); i++) {
        SNP *s = new SNP(*templ.snps[i]);
        this->snps.push_back(s);
    }

    this->comp_cnt  = vector<uint>(templ.comp_cnt);
    this->comp_type = vector<read_type>(templ.comp_type);
    this->loc       = templ.loc;
    this->alleles   = map<string, int>(templ.alleles);
    this->strings   = vector<pair<allele_type, string>>(templ.strings);

}

int
Locus::snp_index(uint col) const
{
    for (uint i = 0; i < this->snps.size(); i++)
        if (this->snps[i]->col == col)
            return i;
    DOES_NOT_HAPPEN;
    return -1;
}

int
Locus::add_consensus(const char *seq)
{
    if (this->con != NULL)
        delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[this->len + 1];
    strcpy(this->con, seq);

    return 0;
}

int
Locus::add_model(const char *seq)
{
    if (this->model != NULL)
        delete [] this->model;

    this->model = new char[this->len + 1];
    strncpy(this->model, seq, this->len);
    this->model[this->len] = '\0';

    return 0;
}

int
Locus::populate_alleles()
{
    vector<SNP *>::iterator  i;
    map<string, int>::iterator j;
    string s;
    uint   k;

    if (this->len > strlen(this->con))
        cerr << "Recorded locus->len: " << this->len << "; consensus length: " << strlen(this->con) << "\n";

    //
    // Is this effective?
    //
    for (uint n = 0; n < this->strings.size(); n++) {
        this->strings[n].first.clear();
        this->strings[n].second.clear();
    }
    this->strings.clear();

    if (this->snps.size() == 0) {
        this->strings.push_back(make_pair("consensus", this->con));
        return 0;
    }

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
        s = this->con;
        k = 0;

        for (i = this->snps.begin(); i != this->snps.end(); i++) {
            if ((*i)->type == snp_type_het && (*i)->col < this->len && k < j->first.length()) {
                s.replace((*i)->col, 1, 1, j->first[k]);
                k++;
            }
        }

        this->strings.push_back(make_pair(j->first, s));
    }

    return 0;
}

int
adjust_snps_for_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, stop, offset, snp_index;

    bp        = 0;
    offset    = 0;
    snp_index = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            offset += dist;
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist;
            while (bp < stop && snp_index < loc->snps.size()) {
                if (loc->snps[snp_index]->col == bp) {
                    loc->snps[snp_index]->col += offset;
                    snp_index++;
                }
                bp++;
            }
            break;
        default:
            break;
        }
    }

    return 0;
}

int
adjust_and_add_snps_for_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, new_bp, stop, snp_cnt;
    SNP   *s;

    bp      = 0;
    new_bp  = 0;
    snp_cnt = loc->snps.size();

    vector<SNP *> snps;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            stop = new_bp + dist;
            while (new_bp < stop) {
                s = new SNP;
                s->col    = new_bp;
                s->type   = snp_type_unk;
                s->rank_1 = 'N';
                snps.push_back(s);
                new_bp++;
            }
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist > snp_cnt ? snp_cnt : bp + dist;
            while (bp < stop) {
                loc->snps[bp]->col = new_bp;
                snps.push_back(loc->snps[bp]);
                bp++;
                new_bp++;
            }
            break;
        default:
            break;
        }
    }

    loc->snps.clear();

    for (uint i = 0; i < snps.size(); i++)
        loc->snps.push_back(snps[i]);

    return 0;
}

int
remove_snps_from_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, new_bp, stop, snp_cnt;

    bp      = 0;
    new_bp  = 0;
    snp_cnt = loc->snps.size();

    vector<SNP *> snps;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            stop = bp + dist;
            while (bp < stop) {
                delete loc->snps[bp];
                bp++;
            }
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist > snp_cnt ? snp_cnt : bp + dist;
            while (bp < stop) {
                loc->snps[bp]->col = new_bp;
                snps.push_back(loc->snps[bp]);
                bp++;
                new_bp++;
            }
            break;
        default:
            break;
        }
    }

    loc->snps.clear();

    for (uint i = 0; i < snps.size(); i++)
        loc->snps.push_back(snps[i]);

    return 0;
}

QLocus::QLocus(const QLocus &other): Locus(other)
{
    for (auto it = other.matches.begin(); it != other.matches.end(); it++) {
        Match *m = *it;
        this->matches.push_back(new Match(*m));
    }
}

QLocus::~QLocus()
{
    for (auto it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance, string cigar)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;
    m->cigar      = cigar;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = "";
    m->dist       = 0;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::clear_matches()
{
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
    this->matches.clear();

    return 0;
}

void
CLocAlnSet::clear() {
    id_= -1;
    aln_pos_.clear();
    ref_ = DNASeq4();
    mpopi_ = NULL;
    reads_.clear();
    reads_per_sample_.clear();
}

size_t CLocReadSet::n_samples() const {
    vector<size_t> n_reads_per_sample (mpopi_->samples().size());
    for (auto& r : reads_)
        ++n_reads_per_sample[r.sample];
    size_t n_samples = 0;
    for (size_t n_s_reads : n_reads_per_sample)
        if (n_s_reads > 0)
            ++n_samples;
    return n_samples;
}

void CLocAlnSet::reinit(int id, const PhyLoc& aln_pos, const MetaPopInfo* mpopi) {
    clear();
    id_ = id;
    aln_pos_ = aln_pos;
    mpopi_ = mpopi;
    reads_per_sample_.resize(mpopi->samples().size());
}

void CLocAlnSet::recompute_consensus() {
    Counts<Nt4> cnts;
    size_t i = 0;
    for (CLocAlnSet::site_iterator site (*this); bool(site); ++site, ++i) {
        site.counts(cnts);
        pair<size_t,Nt4> best_nt = cnts.sorted()[0];
        if (best_nt.first > 0)
            ref_.set(i, best_nt.second);
        else
            ref_.set(i, Nt4::n);
    }
    assert(i == ref_.length());
}

void CLocAlnSet::hard_clip_right_Ns() {
    assert(!ref_.empty());
    //assert(ref_[0] != Nt4::n); // consensus must have been computed, and the first nt is in the cutsite.
    if (ref_[0] == Nt4::n)
        cerr << "WARNING: Some loci seemingly do not start with a cutsite; check your input reads. ("
             << id_ << ")\n";

    size_t to_clip = 0;
    DNASeq4::iterator nt = ref_.end();
    while(nt != ref_.begin() && *--nt == Nt4::n)
        ++to_clip;
    assert(to_clip < ref_.length());

    ref_.resize(ref_.length() - to_clip);

    for (SAlnRead& r : reads_) {
        if (r.cigar.empty() || r.cigar.back().first != 'M' || r.cigar.back().second <= to_clip)
            DOES_NOT_HAPPEN;
        r.seq.resize(ref_.length());
        r.cigar.back().second -= to_clip;
    }
}

void
CLocAlnSet::sort_by_read_name()
{
    sort(this->reads_.begin(), this->reads_.end(),
        [] (const SAlnRead& r1, const SAlnRead& r2) { return r1.name < r2.name; }
        );

    for (auto& s : reads_per_sample_)
        s.clear();
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);
}

void
CLocAlnSet::sort_by_alignment_offset()
{
    auto aln_offset = [] (const Cigar& c) ->size_t {
        assert(!c.empty());
        assert(cigar_is_MDI(c));
        if (c[0].first == 'D') {
            return c[0].second;
        } else if (c[0].first == 'I'
            && c.size() >= 2
            && c[1].first == 'D')
        {
            return c[1].second;
        } else {
            return 0;
        }
    };

    sort(this->reads_.begin(), this->reads_.end(),
        [&aln_offset] (const SAlnRead& r1, const SAlnRead& r2) {
            return aln_offset(r1.cigar) < aln_offset(r2.cigar);
        });

    for (auto& s : reads_per_sample_)
        s.clear();
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);
}

void
CLocAlnSet::merge_paired_reads()
{
    //
    // Sort reads by name. Paired reads should have the same name but end with
    // respectively "/1" and "/2".
    //
    sort_by_read_name();

    // Merge paired reads.
    for (auto r1 = this->reads_.begin(); r1 != this->reads_.end(); ++r1) {
        auto  r2 = r1;
        ++r2;

        if (r2 == this->reads_.end())
            break;

        const string& n1 = r1->name;
        const string& n2 = r2->name;
        const size_t   l = n1.length();

        if (n2.length() == l && l >= 2 &&
            n1[l-2] == '/' && n1[l-1] == '1' &&
            n2[l-2] == '/' && n2[l-1] == '2' &&
            n1.compare(0, l-2, n2, 0, l-2) == 0) {

            // r1 and r2 are paired, merge them.
            //assert(r1->sample == r2->sample); // xxx Fix process_radtags.
            if (r1->sample == r2->sample)
                *r1 = SAlnRead(AlnRead::merger_of(move(*r1), move(*r2)), r1->sample);

            assert(cigar_length_ref(r1->cigar) == ref_.length());

            // Mark r2 for removal and skip it.
            r2->seq.clear();
            ++r1;
            if (r1 == this->reads_.end())
                break;
        }
    }

    // Remove emptied reads.
    reads_.erase(std::remove_if(reads_.begin(),
                                reads_.end(),
                                [](const Read& r) { return r.seq.empty(); }
                                ),
                 reads_.end());

    // Refresh `reads_per_sample_`.
    reads_per_sample_ = vector<vector<size_t>>(mpopi().samples().size());
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);
}

ostream& operator<< (ostream& os, const CLocAlnSet& loc) {
    os << "ref\t.\t" << loc.ref().str();
    for (auto& r : loc.reads_)
        os << "\n" << r.name << "\t" << loc.mpopi().samples()[r.sample].name << "\t" << r.aln();
    return os;
}

void
CLocAlnSet::remove_unmerged_reads(ostream* log)
{
    if (log != NULL)
        *log << "BEGIN unpaired_reads\n";

    for (SAlnRead& r : reads_) {
        if (r.name.back() != 'm') {
            r.seq.clear();
            if (log != NULL)
                *log << "rm_unpaired\t" << r.name << '\n';
        }
    }

    // Remove emptied reads.
    stacks_erase_if(reads_, [](const Read& r){return r.seq.empty();} );

    // Refresh `reads_per_sample_`.
    reads_per_sample_ = vector<vector<size_t>>(mpopi().samples().size());
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);

    if (log != NULL)
        *log << "END unpaired_reads\n";
}

void
CLocAlnSet::remove_pcr_duplicates(vector<size_t>* clone_size_distrib, ostream* log)
{
    if (log != NULL)
        *log << "BEGIN pcr_duplicates\n";

    //
    // Sort reads by (I) sample; (II) insert size; (III, for stability) name.
    //
    for (SAlnRead& r : this->reads_) {
        assert(!r.cigar.empty());
        assert(cigar_is_MDI(r.cigar));
        // assert(r.cigar.front().first == 'M'); // This is the cutsite. //FIXME:
        cigar_canonicalize_MDI_order(r.cigar);
    }
    auto compute_insert_length = [] (const Cigar& c) -> size_t {
        // rem. We assume the asserts above.
        size_t len = cigar_length_ref(c);
        auto last = c.rbegin();
        if (last->first == 'I') {
            len += last->second;
            ++last;
        }
        if (last != c.rend() && last->first == 'D')
            len -= last->second;
        return len;
    };
    sort(this->reads_.begin(), this->reads_.end(),
        [&compute_insert_length](const SAlnRead& r1, const SAlnRead& r2) {
            if (r1.sample != r2.sample)
                return r1.sample < r2.sample;
            size_t i1 = compute_insert_length(r1.cigar);
            size_t i2 = compute_insert_length(r2.cigar);
            if (i1 != i2)
                return i1 < i2;
            return r1.name.compare(r2.name) < 0;
        });

    // Remove reads that have the same insert length.
    vector<SAlnRead>::iterator r = this->reads_.begin();
    assert(r != this->reads_.end());
    size_t len = compute_insert_length(r->cigar);
    while (r != this->reads_.end()) {
        auto group = r;
        size_t group_len = len;
        ++r;
        while (r != this->reads_.end()) {
            len = compute_insert_length(r->cigar);
            if (r->sample != group->sample || len != group_len)
                break;
            ++r;
        }
        size_t clone_size = r - group;
        if (clone_size > 1)
            for (auto r2=group+1; r2!=r; ++r2)
                r2->seq.clear();
        if (clone_size_distrib != NULL) {
            if(clone_size_distrib->size() < clone_size + 1) {
                clone_size_distrib->resize(clone_size + 1);
            }
            ++clone_size_distrib->at(clone_size);
        }
        if (log != NULL) {
            *log << "pcr_dupls\t" << mpopi().samples()[group->sample].name
                    << '\t';
            for (auto r2=group; r2!=r; ++r2)
                *log << ',' << r2->name;;
            *log << '\t' << group_len << '\n';
        }
    }

    // Remove emptied reads.
    stacks_erase_if(reads_, [](const Read& r){return r.seq.empty();} );

    // Refresh `reads_per_sample_`.
    reads_per_sample_ = vector<vector<size_t>>(mpopi().samples().size());
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);

    if (log != NULL)
        *log << "END pcr_duplicates\n";
}

CLocAlnSet
CLocAlnSet::juxtapose(CLocAlnSet&& left, CLocAlnSet&& right, long offset)
{
    assert(left.id() == right.id());
    assert(left.pos() == right.pos());
    assert(&left.mpopi() == &right.mpopi());

    size_t left_ref_len = left.ref().length();
    CLocAlnSet merged (move(left));

    //
    // N.B. Oct. 2017:
    //
    // It actually possible that
    // ``` size_t(-offset) > left.ref().length() ```
    // happens legitimately: if inserts smaller than the read size have been
    // sequenced, the paired-end contig may end upstream of the forward reads,
    // in the barcode/adapter region. In this case, we don't do anything special,
    // as with the current approach the left of the paired-end contig (the
    // `right` matrix) is trimmed anyway.
    //
    // Similarly, it is possible that
    // ``` size_t(-offset) > right.ref().length() ```
    // happens legitimately, for the same reason. In this case, we just move the
    // paired-end reads to the FW object, soft-clipping them entirely.
    //
    // At the moment, we don't do any special trimming of the small-insert reads.
    //
    if (offset < 0 && size_t(-offset) > right.ref().length()) {
        for (SAlnRead& r : right.reads_) {
            r.cigar.assign({{'D', merged.ref().length()}, {'I', r.seq.length()}});
            merged.add(move(r));
        }
        left.clear();
        right.clear();
        return merged;
    }

    // Extend the reference sequence.
    if (offset >= 0) {
        for (long i=0; i<offset; ++i)
            merged.ref_.push_back(Nt4::n);
        merged.ref_.append(right.ref().begin(), right.ref().end());
    } else {
        auto right_itr = right.ref().begin();
        for (long i=0; i<-offset; ++i) {
            assert(right_itr != right.ref().end());
            ++right_itr;
        }
        merged.ref_.append(right_itr, right.ref().end());
    }

    // Extend the left reads.
    for (SAlnRead& r : merged.reads_) {
        cigar_extend_right(r.cigar, right.ref().length() + offset);
        assert(cigar_length_ref(r.cigar) == merged.ref().length());
    }

    // Extend/Trim & add the right reads.
    for (SAlnRead& r : right.reads_) {
        if (offset >= 0) {
           cigar_extend_left(r.cigar, left_ref_len + offset);
        } else {
            cigar_trim_left(r.cigar, -offset);
            assert(cigar_length_query(r.cigar) == r.seq.length());
            cigar_extend_left(r.cigar, left_ref_len);
        }
        merged.add(move(r));
    }
    right.reads_ = vector<SAlnRead>();
    right.reads_per_sample_ = vector<vector<size_t>>(right.mpopi().samples().size());

    left.clear();
    right.clear();
    return merged;
}
