// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __VCF_H__
#define __VCF_H__

#include "constants.h"
#include "nucleotides.h"
#include "utils.h"

/*
 * VcfMeta
 * ==========
 * Represents one line of VCF metainformation.
 */
class VcfMeta {
    string key_;
    string value_;
public:
    VcfMeta(const string& k, const string& p) : key_(k), value_(p) {}

    const string& key() const {return key_;}
    const string& value() const {return value_;}

    struct predefs {
        static const VcfMeta info_AD;
        static const VcfMeta info_AF;
        static const VcfMeta info_DP;
        static const VcfMeta info_NS;

        static const VcfMeta format_AD;
        static const VcfMeta format_DP;
        static const VcfMeta format_GL;
        static const VcfMeta format_GQ;
        static const VcfMeta format_GT;
        static const VcfMeta format_HQ;

        // Custom.
        static const VcfMeta info_loc_strand;
    };
};

/*
 * VcfHeader
 * ==========
 * Stores the contents of a VCF header.
 */
class VcfHeader {
    vector<string> samples_;
    vector<VcfMeta> meta_;
    map<string, size_t> sample_indexes_;

    template<typename OStream> void print(OStream& os) const; // Template because VersatileWriter isn't a proper std::ostream.

public:
    VcfHeader() : samples_(), meta_(), sample_indexes_() {}

    const vector<VcfMeta>& meta() const {return meta_;}
    const vector<string>& samples() const {return samples_;}
    size_t sample_index(const string& sample) const {return sample_indexes_.at(sample);}

    void add_meta(const VcfMeta& m) {meta_.push_back(m);}
    void add_sample(const string& s) {samples_.push_back(s); sample_indexes_.insert({s, samples_.size()-1});}

    // Creates a standard header.
    void add_std_meta(const string& version = "VCFv4.2");

    static const string std_fields;

    friend ostream& operator<< (ostream& os, const VcfHeader& h) {h.print(os); return os;}
    friend VersatileWriter& operator<< (VersatileWriter& w, const VcfHeader& h) {h.print(w); return w;}
};

/*
 * VcfRecord
 * ==========
 * Datastructure to store VCF records
 */
class VcfRecord {
    vector<char> buffer_;
    const char* data() const {return buffer_.data();}

    //size_t chrom_; //0
    size_t pos_;
    //size_t id_;
    size_t allele0_;
    size_t qual_;
    //size_t filters_;
    size_t info0_;
    size_t format0_;
    size_t sample0_;

    void append(const char* str, size_t len) {buffer_.insert(buffer_.end(), str, str+len+1);}
    void append(const string& s) {append(s.c_str(), s.length());}

    bool has_infos() const {return info0_ != SIZE_MAX;}
    bool has_samples() const {return sample0_ != SIZE_MAX;}

    template<typename OStream> void print(OStream& os) const; // Template because VersatileWriter isn't a proper std::ostream.

public:
    VcfRecord()
    : buffer_(), pos_(SIZE_MAX), allele0_(SIZE_MAX), qual_(SIZE_MAX),
      info0_(SIZE_MAX), format0_(SIZE_MAX), sample0_(SIZE_MAX)
    {}
    void assign(const char* rec, size_t len, const VcfHeader& header);

    class iterator {
        const char* p_;
    public:
        iterator(const char* p) : p_(p) {}
        iterator& operator++ () {while(*p_!='\0') ++p_; ++p_; return *this;}
        bool operator!= (iterator other) const {return p_ != other.p_;}
        bool operator== (iterator other) const {return !operator!=(other);}
        const char* operator* () const {return p_;}
    };

    const char* chrom()   const {assert(has_infos()); return data();}
         size_t pos()     const {assert(has_infos()); return atol(data() + pos_) - 1;} // 0-based.
    const char* id()      const {assert(has_infos()); return *++iterator(data() + pos_);}
    const char* allele0() const {assert(has_infos()); return data() + allele0_;}
    const char* qual()    const {assert(has_infos()); return data() + qual_;}
    const char* filters() const {assert(has_infos()); return *++iterator(qual());}
    const char* info0()   const {assert(has_infos()); return data() + info0_;}
    const char* format0() const {assert(has_samples()); return data() + format0_;}
    const char* sample0() const {assert(has_samples()); return data() + sample0_;}

    iterator begin_alleles() const {return iterator(allele0());}
    iterator end_alleles()   const {auto itr=begin_alleles(); ++itr; return **itr=='.' ? itr : iterator(qual());}
    iterator begin_infos()   const {return iterator(info0());}
    iterator end_infos()     const {return *info0()=='.' ? begin_infos() : begin_formats();}
    iterator begin_formats() const {return iterator(format0());}
    iterator end_formats()   const {return *format0()=='.' ? begin_formats() : begin_samples();}
    iterator begin_samples() const {return iterator(sample0());}
    iterator end_samples()   const {return iterator(data() + buffer_.size());}

    size_t count_alleles() const {size_t n=0; for(auto itr=begin_alleles();itr!=end_alleles(); ++itr) ++n; return n;}
    size_t count_infos()   const {size_t n=0; for(auto itr=begin_infos();itr!=end_infos(); ++itr) ++n; return n;}
    size_t count_formats() const {size_t n=0; for(auto itr=begin_formats();itr!=end_formats(); ++itr) ++n; return n;}
    size_t count_samples() const {size_t n=0; for(auto itr=begin_samples();itr!=end_samples(); ++itr) ++n; return n;}

    const char* find_allele(size_t i) const
        {auto a=begin_alleles(); for(size_t j=0; j<i; ++j) {assert(a!=end_alleles()); ++a;} return *a;}
    const char* find_sample(size_t i) const // N.B. Inefficient.
        {auto a=begin_samples(); for(size_t j=0; j<i; ++j) {assert(a!=end_samples()); ++a;} return *a;}

    bool is_monomorphic() const {return **++begin_alleles() == '.';}
    bool is_snp() const;
    size_t index_of_gt_subfield(const char* format_key) const; // SIZE_MAX if not found.

    // Record creation functions.
    void clear();
    void append_chrom(const string& s);
    void append_pos(size_t pos0);
    void append_id(const string& s = ".");
    void append_allele(Nt2 nt);
    void append_allele(const string& s);
    void append_qual(long phred_qual = -1);
    void append_filters(const string& s = ".");
    void append_info(const string& s);
    void append_format(const string& s);
    void append_sample(const string& s);

    friend ostream& operator<< (ostream& os, const VcfRecord& r) {r.print(os); return os;}
    friend VersatileWriter& operator<< (VersatileWriter& w, const VcfRecord& r) {r.print(w); return w;}

    static const char allele_sep = ',';
    static const char filter_sep = ';';
    static const char info_sep = ';';
    static const char format_sep = ':';
    static const char sample_sep = ':';

    struct util {
        static string fmt_info_af(const vector<double>& alt_freqs);
        static string fmt_gt_gl(const vector<Nt2>& alleles, const GtLiks& liks);

        static pair<long,long> parse_gt_gt(const char* gt_str);
        static size_t parse_gt_dp(const char* gt_str, size_t dp_index);
        static Counts<Nt2> parse_gt_ad(const char* gt_str, size_t ad_index, const vector<Nt2>& alleles);
        static uint8_t parse_gt_gq(const char* gt_str, size_t gq_index);
        static GtLiks parse_gt_gl(const char* gt_str, size_t gl_index, const vector<Nt2>& alleles);

        static const char* find_gt_subfield(const char* sample, size_t n);
        static void skip_gt_subfields(const char** start, size_t n);
        static size_t n_possible_genotypes(size_t n_alleles) {return (n_alleles*(n_alleles+1))/2;}
    };
};

/*
 * VcfParser
 * ==========
 */
class VcfParser {
    VersatileLineReader file_;
    VcfHeader header_;

    // Parse the header.
    void read_header();

public:
    VcfParser(const string& path);
    VcfParser(): file_(), header_() {};

    bool next_record(VcfRecord& rec) {
        try {
            const char* line;
            size_t len;
            if (!file_.getline(line, len))
                return false;
            rec.assign(line, len, header_);
            return true;
        } catch (const exception& e) {
            cerr << "Error: At line " << file_.line_number()
                 << " in file '" << file_.path () << "'.\n";
            throw e;
        }
    }

    int open(string &path);
    const VcfHeader& header() const {return header_;};
    const string& path() const {return file_.path();};
    size_t line_number() const {return file_.line_number();}
};

/*
 * VcfWriter
 * ==========
 * (This has become an empty shell...)
 */
class VcfWriter {
private:
    VersatileWriter file_;
    const VcfHeader header_;

public:
    VcfWriter(const string& path, VcfHeader&& header)
        : file_(path), header_(header)
        {file_ << header_;}

    const VcfHeader& header() const {return header_;}
    void write_record(const VcfRecord& r)
        {assert(r.count_samples()==header_.samples().size()); file_ << r;}

    VersatileWriter& file() {return file_;}
};

/*
 * Inline methods.
 * ==========
 */

template<typename OStream>
void VcfHeader::print(OStream& os) const {
    for(const VcfMeta& m : meta())
        os << "##" << m.key() << "=" << m.value() << "\n";

    os << VcfHeader::std_fields;
    if(!samples().empty())
        os << "\tFORMAT";
    for(const string& s : samples())
        os << '\t' << s;
    os << '\n';
}

template<typename OStream>
void VcfRecord::print(OStream& os) const {

    os << chrom()
       << '\t' << pos() + 1
       << '\t' << id();

    // REF & ALT
    iterator itr = begin_alleles();
    iterator end = end_alleles();
    assert(itr!=end);
    os << '\t' << *itr;
    ++itr;
    if (itr == end) {
        assert(strcmp(*itr, ".")==0);
        os << '\t' << '.';
    } else {
        os << '\t' << *itr;
        ++itr;
        for(; itr!=end; ++itr)
            os << allele_sep << *itr;
    }

    //QUAL
    os << '\t' << qual();

    //FILTER
    os << '\t' << filters();

    //INFO
    itr = begin_infos();
    end = end_infos();
    if (itr == end) {
        assert(strcmp(*itr, ".")==0);
        os << '\t' << '.';
    } else {
        os << '\t' << *itr;
        ++itr;
        for(; itr!=end; ++itr)
            os << info_sep << *itr;
    }

    if (begin_samples() != end_samples()) {
        //FORMAT
        itr = begin_formats();
        end = end_formats();
        if (itr == end) {
            assert(strcmp(*itr, ".")==0);
            os << '\t' << '.';
        } else {
            os << '\t' << *itr;
            ++itr;
            for(; itr!=end; ++itr)
                os << format_sep << *itr;
        }

        //SAMPLES
        for (itr=begin_samples(); itr!=end_samples(); ++itr)
            os << '\t' << *itr;
    }

    os << '\n';
}

inline
void VcfRecord::append_chrom(const string& s) {
    assert(buffer_.empty());
    append(s);
}

inline
void VcfRecord::append_pos(size_t pos0) {
    assert(pos_ == SIZE_MAX);
    pos_ = buffer_.size();
    char s[32];
    size_t len = sprintf(s, "%zu", pos0 + 1);
    append(s, len);
}

inline
void VcfRecord::append_id(const string& s) {
    assert(pos_ != SIZE_MAX);
    assert(allele0_ == SIZE_MAX);
    append(s);
}

inline
void VcfRecord::append_allele(Nt2 nt) {
    assert(pos_ != SIZE_MAX);
    assert(buffer_.size() > pos_ + strlen(buffer_.data() + pos_) + 1); // id.
    assert(qual_ == SIZE_MAX);
    if (allele0_ == SIZE_MAX)
        allele0_ = buffer_.size();
    char s[] = {char(nt), '\0'};
    append(s, 1);
}

inline
void VcfRecord::append_allele(const string& s) {
    assert(pos_ != SIZE_MAX);
    assert(buffer_.size() > pos_ + strlen(buffer_.data() + pos_) + 1); // id.
    assert(qual_ == SIZE_MAX);
    if (allele0_ == SIZE_MAX)
        allele0_ = buffer_.size();
    append(s);
}

inline
void VcfRecord::append_qual(long phred_qual) {
    assert(allele0_ != SIZE_MAX);
    assert(qual_ == SIZE_MAX);
    qual_ = buffer_.size();
    if (phred_qual >= 0) {
        char s[32];
        size_t len = sprintf(s, "%ld", phred_qual);
        append(s, len);
    } else {
        append(".", 1);
    }
}

inline
void VcfRecord::append_filters(const string& s) {
    assert(qual_ != SIZE_MAX);
    assert(info0_ == SIZE_MAX);
    append(s);
}

inline
void VcfRecord::append_info(const string& s) {
    assert(qual_ != SIZE_MAX);
    assert(buffer_.size() > qual_ + strlen(buffer_.data() + qual_) + 1); // filters.
    assert(format0_ == SIZE_MAX);
    if (info0_ == SIZE_MAX)
        info0_ = buffer_.size();
    append(s);
}

inline
void VcfRecord::append_format(const string& s) {
    assert(info0_ != SIZE_MAX);
    assert(sample0_ == SIZE_MAX);
    if (format0_ == SIZE_MAX)
        format0_ = buffer_.size();
    append(s);
}

inline
void VcfRecord::append_sample(const string& s) {
    assert(format0_ != SIZE_MAX);
    if (sample0_ == SIZE_MAX)
        sample0_ = buffer_.size();
    append(s);
}

inline
void VcfRecord::clear() {
    buffer_.resize(0);
    pos_ = SIZE_MAX;
    allele0_ = SIZE_MAX;
    qual_ = SIZE_MAX;
    info0_ = SIZE_MAX;
    format0_ = SIZE_MAX;
    sample0_ = SIZE_MAX;
}

inline
size_t VcfRecord::index_of_gt_subfield(const char* key) const {
    size_t i=0;
    for (auto f=begin_formats(); f!=end_formats(); ++f) {
        if (strcmp(*f, key) == 0)
            return i;
        ++i;
    }
    return SIZE_MAX;
}

inline
pair<long,long> VcfRecord::util::parse_gt_gt(const char* sample)
{ try {
    assert(sample != NULL && sample[0] != '\0');
    if (*sample == '.')
        return {-1, -1};
    char* end;
    long first = strtol(sample, &end, 10);
    if (end == sample || (*end != '/' && *end != '|'))
        throw exception();
    sample = end;
    ++sample;
    if (*sample == '.')
        // Incomplete, e.g. "1/."
        return {-1, -1};
    long second = strtol(sample, &end, 10);
    if (end == sample || (*end != ':' && *end != '\0'))
        throw exception();
    if (first < 0 || second < 0)
        throw exception();
    return {first, second};
} catch (exception&) {
    cerr << "Error: Malformed VCF sample field '" << sample << "'.\n";
    throw;
}}

inline
bool VcfRecord::is_snp() const {
    iterator a = begin_alleles();
    const iterator end = end_alleles();
    if (strlen(*a) > 1)
        return false;
    ++a;
    if (a == end)
        return false;
    assert(**a != '.');
    if (**a == '*')
        return false;
    for (; a != end; ++a)
        if (strlen(*a) > 1)
            return false;
    return true;
}

inline
const char* VcfRecord::util::find_gt_subfield(const char* sample, size_t n) {
    const char* subf = sample;
    for(size_t i=0; i<n; ++i) {
        subf = strchr(subf, ':');
        if (subf == NULL)
            break;
        ++subf;
    }
    return subf;
}

#endif // __VCF_H__
