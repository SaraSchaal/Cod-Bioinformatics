#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "constants.h"
#include "DNASeq4.h"
#include "stacks.h"
#include "aln_utils.h"

class Alignment {
    const DNASeq4* seq_;
    const Cigar* cig_;

public:
    Alignment(const DNASeq4& seq, const Cigar& cigar) : seq_(&seq), cig_(&cigar) {
        assert(cigar_length_query(*cig_) == seq_->length());
        assert(check_cigar_ops());
    }

    // N.B. Inefficient; prefer iteration.
    Nt4 operator[] (size_t ref_i) const;

    friend ostream& operator<< (ostream& os, const Alignment& aln)
        {for(iterator it (aln); it; ++it) os << char(*it); return os;}

private:
    bool check_cigar_ops() const;

public:
    // Iterator.
    // We have to use a range-style iterator to skip insertions (as we can't
    // peek at the next CIGAR operation if we don't know if we have reached the
    // end or not / we can't dereference cig_.end()).
    class iterator {
        Cigar::const_iterator cig_it_;
        Cigar::const_iterator cig_past_;
        size_t pos_; // Position in the current cigar op.
        DNASeq4::iterator seq_it_;
        DNASeq4::iterator seq_past_; // For debugging purposes; only used in assert's.

    public:
        iterator(const Alignment& a)
            : cig_it_(a.cig_->begin()), cig_past_(a.cig_->end()), pos_(0), seq_it_(a.seq_->begin()), seq_past_(a.seq_->end())
            {skip_insertion();}
        iterator& operator++ ();
        iterator& operator+= (size_t n);
        operator bool() const {return cig_it_ != cig_past_;}

        Nt4 operator* () const {if (cig_it_->first=='M') return *seq_it_; else {assert(cig_it_->first=='D'); return Nt4::n;}}
        char curr_op() const {return cig_it_->first;} // N.B. This is never 'I' (as it is skipped).

    private:
        void skip_insertion();
    };
};

struct AlnRead : Read {
    Cigar cigar;
    AlnRead() : Read() {}
    AlnRead(Read&& r, Cigar&& c) : Read(move(r)), cigar(move(c)) {}

    Alignment aln() const {return Alignment(seq, cigar);}

    static AlnRead merger_of(const AlnRead& r1, const AlnRead& r2);
};

//
// ==================
// Inline definitions
// ==================
//

inline
Nt4 Alignment::operator[] (size_t ref_i) const
{
    size_t seq_i = 0;
    for (auto op=cig_->begin(); op!=cig_->end(); ++op) {
        if (op->first == 'M') {
            if (ref_i < op->second)
                // This is the relevant cigar operation.
                return (*seq_)[seq_i+ref_i];

            // Consume ref & seq.
            seq_i += op->second;
            ref_i -= op->second;

        } else if (op->first == 'D') {
            if (ref_i < op->second)
                // This is the relevant cigar operation.
                return Nt4::n;

            // Consume ref.
            ref_i -= op->second;
        } else if (op->first == 'I') {
            // Consume seq.
            seq_i += op->second;
        } else {
            DOES_NOT_HAPPEN;
        }
    }
    // `ref_i` wasn't entirely consumed.
    throw std::out_of_range("Alignment::op[]: out_of_range");
    return Nt4::n;
}

inline
Alignment::iterator& Alignment::iterator::operator++ ()
{
    assert(cig_it_ != cig_past_);

    //
    // If the current cigar operation is M, advance in the sequence.
    // (rem. The current operation is never I, as they're skipped.)
    //
    if (cig_it_->first == 'M') {
        assert(seq_it_ != seq_past_);
        ++seq_it_;
    }

    //
    // Advance in the cigar.
    //
    ++pos_;
    if (pos_ == cig_it_->second) {
        // Enter the next CIGAR operation.
        pos_ = 0;
        ++cig_it_;
        skip_insertion();
    }
    // Upon reaching the end of the CIGAR, check that the entire sequence was
    // also consumed.
    assert(cig_it_ == cig_past_ ? !(seq_it_ != seq_past_) : true);

    return *this;
}

inline
Alignment::iterator& Alignment::iterator::operator+= (size_t n) {
    // We could block-increment according to the cigar, but this is not implemented
    // at the moment, i.e. we just increment n times.
    for (size_t i=0; i<n; ++i)
        ++*this;
    return *this;
}

inline
void
Alignment::iterator::skip_insertion()
{
    if (cig_it_ != cig_past_ && cig_it_->first == 'I') {
        // Op is I; skip this insertion.
        for (size_t i=0; i<cig_it_->second; ++i) {
            assert(seq_it_ != seq_past_);
            ++seq_it_;
        }
        ++cig_it_;
    }
}

inline bool
Alignment::check_cigar_ops() const
{
    for (auto& op : *cig_)
        if (op.first != 'M' && op.first != 'D' && op.first != 'I')
            // Oops, the class only knows M, I and D.
            return false;

    return true;
}

inline AlnRead
AlnRead::merger_of(const AlnRead& r1, const AlnRead& r2)
{
    assert(cigar_length_ref(r1.cigar) == cigar_length_ref(r2.cigar));

    // name.
    string name (r1.name);
    if (name[name.length()-2] == '/' && (name.back() == '1' || name.back() == '2'))
        name.back() = 'm';
    else
        name += "/m";

    // sequence & cigar.
    // xxx We use push_back(); it would be faster to append (but it would require a more elaborate handling of cigars).
    DNASeq4 seq;
    Cigar cig;
    auto a1 = Alignment::iterator(r1.aln());
    auto a2 = Alignment::iterator(r2.aln());
    auto increment_cigar_op = [&cig](char op){
        if (cig.empty() || cig.back().first != op)
            cig.push_back({op, 1});
        else
            ++cig.back().second;
    };
    while(a1) {
        assert(a2); // a1 & a2 must become false together.
        char o1 = a1.curr_op();
        char o2 = a2.curr_op();
        assert((o1 == 'M' || o1 == 'D') && (o2 =='M' || o2 == 'D'));
        if (o1 == 'M' && o2 == 'M') {
            Nt4 nt1 = *a1;
            Nt4 nt2 = *a2;
            if (nt1 == nt2)
                seq.push_back(nt1);
            else
                seq.push_back(Nt4::n);
            increment_cigar_op('M');
        } else if (o1 == 'D' && o2 == 'D') {
            increment_cigar_op('D');
        } else if (o1 == 'M' && o2 == 'D') {
            Nt4 nt1 = *a1;
            seq.push_back(nt1);
            increment_cigar_op('M');
        } else if (o1 == 'D' && o2 == 'M') {
            Nt4 nt2 = *a2;
            seq.push_back(nt2);
            increment_cigar_op('M');
        }
        ++a1;
        ++a2;
    }
    assert(!a2);

    return AlnRead(Read(move(seq), move(name)), move(cig));
}

#endif
