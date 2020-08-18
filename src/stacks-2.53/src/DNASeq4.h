#ifndef DNASEQ4_H
#define DNASEQ4_H

#include "constants.h"
#include "nucleotides.h"

// A sequence of nucleotides on a uint64_t, where the first nucleotide uses the
// low bits.
template<class Nt>
class NtArray {
    uint64_t a_;

public:
    NtArray() : a_(0) {}
    NtArray(uint64_t x) : a_(x) {}
    NtArray(const NtArray<Nt>& other) : a_(other.a_) {}
    NtArray<Nt>& operator= (const NtArray<Nt>& other) {a_ = other.a_; return *this;}

    void set(size_t i, Nt nt) {a_ |= uint64_t(size_t(nt)) << (i*Nt::nbits);}
    void clear(size_t i) {a_ &= ~(lowbits << (i*Nt::nbits));}

    Nt operator[] (size_t i) const {return Nt(size_t((a_ >> (i*Nt::nbits)) & lowbits));}
    bool operator== (const NtArray<Nt>& other) const {return a_ == other.a_;}
    bool operator<  (const NtArray<Nt>& other) const {return a_ < other.a_;}
    friend struct std::hash<NtArray>;

    // Methods for kmers
    void push_front(Nt nt) {a_ <<= Nt::nbits; a_ |= uint64_t(size_t(nt));}
    void pop_front() {a_ >>= Nt::nbits;}

    static const size_t n_nts = sizeof a_ * 8 / Nt::nbits;

private:
    static const uint64_t lowbits = (1 << Nt::nbits) - 1;
};

// A dinucleotide coded on one byte, where the first nucleotide uses the high bits.
// This is compatible with BAM/HTSLIB.
class DiNuc {
    uchar x_;

public:
    DiNuc() : x_(0) {}
    DiNuc(const DiNuc& other) : x_(other.x_) {}
    DiNuc(uchar x) : x_(x) {}
    DiNuc(Nt4 nt1, Nt4 nt2) : x_(0) {set_first(nt1); set_second(nt2);}
    DiNuc(char nt1, char nt2) : DiNuc(Nt4(nt1), Nt4(nt2)) {}
    DiNuc& operator= (const DiNuc& other) {x_ = other.x_; return *this;}

    Nt4 first() const {return Nt4(size_t(x_ >>4));}
    Nt4 second() const {return Nt4(size_t(x_ & 15));}
    void first(Nt4 nt) {clear_first(); set_first(nt);}
    void second(Nt4 nt) {clear_second(); set_second(nt);}

    bool operator== (const DiNuc& other) const {return x_ == other.x_;}
    bool operator<  (const DiNuc& other) const {return x_ < other.x_;}

private:
    void set_first(Nt4 nt) {x_ |= uchar(size_t(nt)) <<4;}
    void set_second(Nt4 nt) {x_ |= uchar(size_t(nt));}

    void clear_first() {x_ &= 15;}
    void clear_second() {x_ &= ~15;}

public:
    uchar x() const {return x_;} // c.f. hash<DNASeq4>
};

// A sequence of 4-bits A, C, G, T and N.
// Uses DiNuc and is thus compatible with BAM/HTSLIB.
class DNASeq4 {
    size_t l_;
    vector<DiNuc> v_;

public:
    class iterator;

    DNASeq4() : l_(0), v_() {}
    explicit DNASeq4(const DNASeq4& other) : l_(other.l_), v_(other.v_) {}
    DNASeq4(DNASeq4&& other) noexcept : l_(other.l_), v_(move(other.v_)) {other.clear();}
    DNASeq4(const char* s, size_t len);
    DNASeq4(const uchar* arr, size_t len) : l_(len), v_(arr, arr+len/2+len%2) {} // `len` is the length of the sequence (not of the array)
    DNASeq4(const string& s) : DNASeq4(s.c_str(), s.size()) {}
    DNASeq4(size_t len) : l_(len), v_(l_/2+l_%2, DiNuc(Nt4::n, Nt4::n)) {if(l_%2) v_.back().second(Nt4::$);}
    DNASeq4& operator= (const DNASeq4& other) {l_ = other.l_; v_ = other.v_; return *this;}
    DNASeq4& operator= (DNASeq4&& other) {l_ = other.l_; v_ = move(other.v_); other.clear(); return *this;}

    size_t length() const {return l_;}
    bool empty() const {return l_ == 0;}
    void set(size_t i, Nt4 nt) {i%2==0 ? v_[i/2].first(nt) : v_[i/2].second(nt);}
    void clear() {l_ = 0; v_.clear();}
    void resize(size_t len);
    void reserve(size_t len) {v_.reserve(len/2+len%2);}
    void push_back(Nt4 nt) {++l_; if (l_%2) v_.push_back(DiNuc(nt,Nt4::$)); else v_.back().second(nt);}
    void append(iterator first, iterator past);

    DNASeq4 rev_compl() const;
    string str() const;

    void remove_Ns();
    void shift_Ns_towards_the_end();

    Nt4 operator[] (size_t i) const {return i%2==0 ? v_[i/2].first() : v_[i/2].second();}
    bool  operator== (const DNASeq4& other) const {return l_ == other.l_ && v_ == other.v_;}
    bool  operator<  (const DNASeq4& other) const {return l_ < other.l_ ? true : v_ < other.v_;}
    friend struct std::hash<DNASeq4>;
    friend ostream& operator<< (ostream& os, const DNASeq4& seq) {for (Nt4 nt : seq) os << nt; return os;}

    // Iterator.
    class iterator {
        vector<DiNuc>::const_iterator vi_;
        bool first_;

    public:
        iterator(vector<DiNuc>::const_iterator vi, bool f) : vi_(vi), first_(f) {}
        bool operator!= (iterator other) const {return ! (vi_ == other.vi_? first_ == other.first_ : false);}
        bool operator== (iterator other) const {return !operator!=(other);}
        iterator& operator++ () {if (first_) {first_ = false;} else {++vi_; first_ = true;} return *this; }
        iterator& operator-- () {if (first_) {--vi_; first_ = false;} else {first_ = true;} return *this; }
        size_t operator- (iterator other) {return 2*(vi_ - other.vi_) + size_t(other.first_) - size_t(first_);}

        // Get the (Nt4) nucleotide.
        Nt4 operator* () const {return first_ ? vi_->first() : vi_->second();}
    };
    iterator begin() const {return iterator(v_.begin(), true);}
    iterator end()   const {return length()%2==0 ? iterator(v_.end(), true) : iterator(--v_.end(), false);}

    // Methods to allow to memcpy into a htslib bam1_t.
    size_t nbytes() const {return v_.size() * sizeof (DiNuc);}
    const uchar* vdata() const {return (uchar*) v_.data();}
};

//
// Inline definitions.
// ==========
//

namespace std { template<class Nt>
struct hash<NtArray<Nt>> { size_t operator() (const NtArray<Nt>& a) const {
    return hash<uint64_t>()(a.a_);
}};}

namespace std {
template<>
struct hash<DNASeq4> {
    size_t operator() (const DNASeq4& s) const {
        size_t x = static_cast<size_t>(14695981039346656037ULL);
        for (const DiNuc& d : s.v_) {
            x ^= static_cast<size_t>(d.x());
            x *= static_cast<size_t>(1099511628211ULL);
        }
        return x;
    }
};
}

#endif //DNASEQ4_h
