// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __BAMI_H__
#define __BAMI_H__

//
// Code to parse binary BAM format. This format is created for
// reads that have been aligned to a reference genome.
//

#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include "htslib/sam.h"

#include "stacks.h"
#include "input.h"
#include "DNASeq4.h"

typedef vector<pair<char,uint>> Cigar;

// Trivial hash table giving the binary value (0-9) of a cigar character.
// htslib doesn't define this (c.f. the use of `bam_hdr_t::cigar_tab` in sam.c).
// However, it does define the reverse variable in sam.h:
// `#define BAM_CIGAR_STR "MIDNSHP=XB"`
extern const uint32_t cigar_c2i[128];

int bam_find_start_bp(int, strand_type, const Cigar&);
int bam_edit_gaps(const Cigar&, char*);

class BamRecord {
    bam1_t* r_;

    BamRecord& operator= (const BamRecord&) = delete;

public:
    BamRecord()                              : r_(NULL) {}
    BamRecord(BamRecord&& other) noexcept    : BamRecord() {std::swap(r_, other.r_);} // noexcept necessary for std::vector to use it...
    explicit BamRecord(const BamRecord&)     : BamRecord() {DOES_NOT_HAPPEN;} // Necessary for std::vector to compile...
    ~BamRecord()                             {destroy();}
    BamRecord& operator= (BamRecord&& other) {std::swap(r_, other.r_); return *this;}
    void reinit()                            {destroy(); r_=bam_init1(); r_->data=NULL; r_->m_data=0;}
    void destroy()                           {if (r_!=NULL) {bam_destroy1(r_); r_=NULL;}}
    bool empty() const                       {return r_==NULL;}

    // Access the fields.
    const char* qname()     const {return bam_get_qname(r_);}
    uint8_t     qname_len() const {return r_->core.l_qname-1;}
    uint16_t    flag()      const {return r_->core.flag;}
    int32_t     chrom()     const {return r_->core.tid;}
    int32_t     pos()       const {return r_->core.pos;}
    uint8_t     mapq()      const {return r_->core.qual;}
    Cigar       cigar()     const;
    // (rnext, pnext, tlen)
    DNASeq4     seq()       const {return DNASeq4((uchar*) bam_get_seq(r_), r_->core.l_qseq);}
    // (qual, aux)

    // Flags.
    bool is_unmapped()      const {return r_->core.flag & BAM_FUNMAP;}
    bool is_secondary()     const {return r_->core.flag & BAM_FSECONDARY;}
    bool is_supplementary() const {return r_->core.flag & BAM_FSUPPLEMENTARY;}
    bool is_primary()       const {return ! (r_->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY));}
    bool is_rev_compl()     const {return r_->core.flag & BAM_FREVERSE;}
    bool is_read1()         const {return r_->core.flag & BAM_FREAD1;}
    bool is_read2()         const {return r_->core.flag & BAM_FREAD2;}
    // (etc.)

    void assign(
            const string& name,
            uint16_t flg,
            int32_t chr_index,
            int32_t aln_pos,
            const Cigar& cig,
            const DNASeq4& seq,
            size_t read_group
            );

    // Access to the underlying hts record.
    const bam1_t*   hts()         const {return r_;}
    bam1_t*         hts()               {return r_;}
    int32_t         hts_l_seq()   const {return r_->core.l_qseq;}
    const uchar*    hts_seq()     const {return bam_get_seq(r_);}
    uint32_t        hts_n_cigar() const {return r_->core.n_cigar;}
    const uint32_t* hts_cigar()   const {return bam_get_cigar(r_);}

    // Other convenience functions.
    static uint32_t cig_op_t(uint32_t op) {return op & BAM_CIGAR_MASK;}
    static uint32_t cig_op_len(uint32_t op) {return op >> BAM_CIGAR_SHIFT;}
    static uint32_t cig_ref_len(const uint32_t* cigar, uint32_t n_cigar);
    static uint32_t cig_ref_len(const BamRecord& r) {return cig_ref_len(r.hts_cigar(), r.hts_n_cigar());}
    const char* read_group() const;

private:
    // Moves the given pointer to the start of the next AUX field.
    static void skip_one_aux(const uint8_t** ptr);
};

class BamHeader {
    bam_hdr_t* h_;
    void init_sdict(); // (So that BamHeader::chrom_index can be const.)

    BamHeader& operator= (BamHeader&) = delete;

public:
    BamHeader()                             : h_(NULL) {}
    BamHeader(const char* text, size_t len);
    BamHeader(const string& text)           : BamHeader(text.c_str(), text.length()) {}
    explicit BamHeader(const BamHeader& other) : BamHeader(other.text(), strlen(other.text())) {}
    BamHeader(BamHeader&& other)            : BamHeader() {std::swap(h_, other.h_);}
    ~BamHeader()                            {destroy();}
    BamHeader& operator=(BamHeader&& other) {std::swap(h_, other.h_); return *this;}
    void reinit()                           {destroy(); h_ = bam_hdr_init();}
    void reinit(htsFile* bam_f)             {destroy(); h_ = bam_hdr_read(bam_f->fp.bgzf); init_sdict();}
    void destroy()                          {if (h_!=NULL) {bam_hdr_destroy(h_); h_=NULL;}}
    bool empty() const                      {return h_==NULL;}

    const char* text()                   const {assert(h_); return h_->text;}
    size_t      n_ref_chroms()           const {assert(h_); return h_->n_targets;}
    const char* chrom_str(int32_t index) const;
    size_t      chrom_len(int32_t index) const;
    int32_t     chrom_index(const char* chrom_str) const;

    const bam_hdr_t* hts() const {return h_;}
          bam_hdr_t* hts()       {return h_;}

    static void check_same_ref_chroms(const BamHeader& h1, const BamHeader& h2);
};

class Bam: public Input {
    htsFile*  bam_fh;
    BamHeader hdr;
    bool eof_;

    size_t  n_records_read_;
    int32_t prev_chrom_;
    int32_t prev_pos_;

    Bam(Bam&&) = delete; // (This is just not implemented at the moment.)

public:
    Bam(const char* path);
    Bam(const string& path) : Bam(path.c_str()) {}
    Bam(const string& path, BamHeader&& header);
    ~Bam() {if (bam_fh) hts_close(bam_fh);}

    const BamHeader& h() const {return hdr;}

    bool next_record(BamRecord& rec);
    bool next_record_ordered(BamRecord& rec);
    bool eof() const {return eof_;}
    void close() { assert(bam_fh); if(hts_close(bam_fh) != 0) throw ios::failure("hts_close"); bam_fh = NULL;}
    size_t n_records_read() const {return n_records_read_;}

    Seq *next_seq();
    int  next_seq(Seq&);

    void write(const BamRecord& rec) {if(bam_write1(bam_fh->fp.bgzf, rec.hts()) < 0) throw ios::failure("bam_write1");}

private:
    static void check_open(const htsFile* bam_f, const string& path);
};

//
// Inline definitions
// ----------
//

inline
bool Bam::next_record(BamRecord& rec) {
    if (rec.empty())
        rec.reinit();
    int rv = bam_read1(bam_fh->fp.bgzf, rec.hts());
    if (rv == -1) {
        eof_ = true;
        rec.destroy();
        return false;
    } else if (rv < -1) {
        cerr << "Error: while reading BAM file '" << this->path << "'.\n";
        throw ios::failure("hts::bam_read1");
    }
    ++n_records_read_;
    return true;
}

inline
bool Bam::next_record_ordered(BamRecord& rec) {
    if (!next_record(rec))
        return false;
    // Check the ordering.
    if (rec.chrom() >= 0) {
        if (prev_chrom_ >= 0) {
            if (rec.chrom() < prev_chrom_
                || (rec.pos() < prev_pos_ && rec.chrom() == prev_chrom_)
            ){
                cerr << "Error: BAM file is not properly sorted; " << n_records_read_
                     << "th record '" << rec.qname() << "' at "
                     << hdr.chrom_str(rec.chrom()) << ':' << rec.pos()
                     << " should come before previously seen position " << hdr.chrom_str(prev_chrom_)
                     << ':' << prev_pos_ << ".\n";
                throw exception();
            }
        } else {
            if (rec.chrom() >= 0) {
                cerr << "Error: BAM file is not properly sorted; " << n_records_read_
                     << "th record '" << rec.qname() << "' at "
                     << hdr.chrom_str(rec.chrom()) << ':' << rec.pos()
                     << "' should come before unmapped records.\n";
                throw exception();
            }
        }
    }
    prev_chrom_ = rec.chrom();
    prev_pos_ = rec.pos();
    return true;
}

inline
Cigar BamRecord::cigar() const {
    Cigar cig;
    for (size_t i = 0; i < hts_n_cigar(); ++i) {
        uint32_t op = hts_cigar()[i];
        char op_t;
        switch(cig_op_t(op)) {
        case BAM_CMATCH:
            op_t = 'M';
            break;
        case BAM_CEQUAL:
            op_t = '=';
            break;
        case BAM_CDIFF:
            op_t = 'X';
            break;
        case BAM_CINS:
            op_t = 'I';
            break;
        case BAM_CDEL:
            op_t = 'D';
            break;
        case BAM_CSOFT_CLIP:
            op_t = 'S';
            break;
        case BAM_CREF_SKIP:
            op_t = 'N';
            break;
        case BAM_CHARD_CLIP:
            op_t = 'H';
            break;
        case BAM_CPAD:
            op_t = 'P';
            break;
        default:
            cerr << "Warning: Unknown CIGAR operation at record '" << bam_get_qname(r_)
                 << "(binary value " << cig_op_t(op) << ").\n";
            op_t = '?';
            break;
        }
        cig.push_back({op_t, cig_op_len(op)});
    }
    return cig;
}

inline
uint32_t BamRecord::cig_ref_len(const uint32_t* cigar, uint32_t n_cigar) {
    int32_t ref_len = 0;
    for (size_t i = 0; i < n_cigar; ++i) {
        uint32_t op = cigar[i];
        switch(cig_op_t(op)) {
        case BAM_CINS:
        case BAM_CSOFT_CLIP:
            break;
        default:
            ref_len += cig_op_len(op);
            break;
        }
    }
    return ref_len;
}

inline
const char* BamRecord::read_group() const {
    const uint8_t* ptr = bam_get_aux(r_);
    while (ptr < bam_get_aux(r_) + bam_get_l_aux(r_) - 3)
        if (strncmp((const char*)ptr, "RGZ", 3) == 0)
            return (const char*)ptr+3;
        else
            skip_one_aux(&ptr);

    return NULL;
}

inline
void BamRecord::skip_one_aux(const uint8_t** ptr) {
    using namespace std;

    // The library doesn't provide much for handling the AUX fields.

    // Make sure that the tag matches [A-Za-z][A-Za-z0-9].
    if (!isalpha(**ptr) || !isalnum(*(*ptr+1))) {
        cerr << "Warning: Illegal BAM AUX tag '"
             << *(const char*)*ptr << *(const char*)(*ptr+1) << "' (ASCII "
             << **ptr << "," << *(*ptr+1) << ").\n" << flush;
    }
    *ptr += 2;

    uint8_t t = **ptr; // byte 3 of the field gives the type
    *ptr += 1;

    if (t == 'i'        // int32_t
            || t == 'I' // uint32_t
            || t == 'f' // float
            ) {
        *ptr += 4;
    } else if (t == 'Z' // string
            || t == 'H' // hex string
            ) {
        *ptr += strlen((const char*)*ptr) + 1;
    } else if (t == 'A' // character
            || t == 'c' // int8_t
            || t == 'C' // uint8_t
            ) {
        *ptr += 1;
    } else if (t == 's' // int16_t
            || t == 'S' // uint16_t
            ) {
        *ptr += 2;
    } else if (t == 'B') { // array
        *ptr += 1;
        uint8_t t2 = **ptr; // byte 4 gives the array contents type
        *ptr += 1;
        int32_t len = *(int32_t*)*ptr; // bytes 5-8 give the array length
        *ptr += 4;

        // Type should be an integer type or float.
        if (t2 == 'c' || t2 == 'C') {
            *ptr += len;
        } else if (t2 == 's' || t2 == 'S') {
            *ptr += 2 * len;
        } else if (t2 == 'i' || t2 == 'I' || t2 == 'f') {
            *ptr += 4 * len;
        } else {
            cerr << "Error: Unexpected BAM AUX field array type '" << (char)t2 << "' (ASCII " << t2 << ").\n";
            throw exception();
        }
    } else {
        cerr << "Error: Unexpected BAM AUX field type '" << (char)t << "' (ASCII " << t << ").\n";
        throw exception();
    }
}

inline
const char* BamHeader::chrom_str(int32_t index) const
{
    assert(h_);
    if (index < 0 || index >= h_->n_targets)
        throw std::out_of_range("BamHeader::chrom_str");
    return h_->target_name[index];
}

inline
size_t BamHeader::chrom_len(int32_t index) const
{
    assert(h_);
    if (index < 0 || index >= h_->n_targets)
        throw std::out_of_range("BamHeader::chrom_len");
    return h_->target_len[index];
}

inline
int32_t BamHeader::chrom_index(const char* chrom_str) const
{
    assert(h_);
    assert(h_->sdict);
    int32_t index = bam_name2id(h_, chrom_str);
    if (index < 0)
        throw std::out_of_range("BamHeader::chrom_index");
    return index;
}

#endif // __BAMI_H__
