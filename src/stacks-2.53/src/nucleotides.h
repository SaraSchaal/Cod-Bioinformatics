#ifndef NUCLEOTIDES_H
#define NUCLEOTIDES_H

#include "constants.h"

class Nt2;

//
// Nt4: A nucleotide coded on 4 bits.
// These definitions are compatible with those of htslib--htslib supports
// partially ambiguous nucleotides ('R', etc.) but we convert everything to 15
// (i.e. 0xF, 'N').
//
class Nt4 {
    size_t nt_;

public:
    Nt4() : nt_(Nt4::n) {}
    Nt4(char c) : nt_(from_ch[size_t(c)]) {}
    Nt4(size_t i) : nt_(i) {}
    Nt4(const Nt4& other) : nt_(other.nt_) {}
    Nt4(const Nt2 nt2);
    Nt4& operator=(const Nt4& other) {nt_ = other.nt_; return *this;}

    Nt4 rev_compl() const {return rev_compl_[size_t(nt_)];}
    bool is_acgt() const {return *this == a || *this == c || *this == g || *this == t;}
    bool is_acgtn() const {return is_acgt() || *this == n ;}
    size_t index() const {return to_index[size_t(nt_)];}
    explicit operator size_t () const {return nt_;}
    explicit operator int () const {return nt_;}
    explicit operator char () const {return to_ch[nt_];}
    bool operator== (Nt4 other) const {return nt_ == other.nt_;}
    bool operator!= (Nt4 other) const {return !operator==(other);}
    bool operator< (Nt4 other) const {return nt_ < other.nt_;}
    friend ostream& operator<< (ostream& os, Nt4 nt) {os << char(nt); return os;}

    static constexpr size_t nbits = 4;
    static const Nt4 $; // 0000 (0)
    static const Nt4 a; // 0001 (1)
    static const Nt4 c; // 0010 (2)
    static const Nt4 g; // 0100 (4)
    static const Nt4 t; // 1000 (8)
    static const Nt4 n; // 1111 (15)
    static const array<Nt4,5> all; // All of the above.

    static constexpr size_t max() {return (1 << nbits) - 1;}

private:
    Nt4(int i) : nt_(i) {}

    // Trivial ASCII-like hash table giving the 4-bits value of a nucleotide letter.
    // e.g. ch_to_nt [ (int)'G' ] == 4
    // Adapted from `htslib::seq_nt16_table` (hts.cc).
    static const Nt4 from_ch[256];

    // Trivial hash table to convert Nt2 nucleotides to Nt4.
    static const Nt4 from_nt2[4];

    // Trivial hash table giving the nucleotide letter of a 4-bits value.
    // e.g. nt_to_ch[4] == 'C'
    static const char to_ch[16];

    // Table giving the reverse complement of the nucleotide.
    static const Nt4 rev_compl_[16];

    // Convert $,a,c,g,t,n into a 0-5 index for arrays of objects indexed by possible Nt4 characters.
    static const size_t to_index[16];
};

//
// Nt2: A nucleotide coded on 2 bits.
//
class Nt2 {
    size_t nt_;

public:
    Nt2() : nt_(Nt2::a) {}
    Nt2(char c) : nt_(from_ch[size_t(c)]) {assert(Nt4(c).is_acgt());}
    Nt2(Nt4 nt4) : nt_(from_nt4[size_t(nt4)]) {assert(nt4.is_acgt());}
    Nt2(size_t i) : nt_(i) {}
    Nt2(int i) : nt_(i) {}
    Nt2(const Nt2& other) : nt_(other.nt_) {}
    Nt2& operator=(const Nt2& other) {nt_ = other.nt_; return *this;}

    Nt2 rev_compl() const {return rev_compl_[nt_];}

    explicit operator size_t () const {return nt_;}
    explicit operator int () const {return nt_;}
    explicit operator char () const {return to_ch[nt_];}
    bool operator== (Nt2 other) const {return nt_ == other.nt_;}
    bool operator!= (Nt2 other) const {return !operator==(other);}
    bool operator< (Nt2 other) const {return nt_ < other.nt_;}
    friend ostream& operator<< (ostream& os, Nt2 nt) {os << char(nt); return os;}

    static constexpr size_t nbits = 2;
    static const Nt2 a;
    static const Nt2 c;
    static const Nt2 g;
    static const Nt2 t;
    static const array<Nt2,4> all;

    static constexpr size_t max() {return (1 << nbits) - 1;}

private:
    static const Nt2 from_ch[256];
    static const Nt2 from_nt4[16];
    static const char to_ch[4];

    static const Nt2 rev_compl_[4];
};

//
// Counts: A class to store nucleotide counts.
// e.g. Counts<Nt2> is a std::array of {n_A, n_C, n_G, n_T}.
//
template<typename Nt>
class Counts {
    // Array of counts, containing the count of A's at index Nt::a,of C's at
    // index Nt::c, etc.
    array<size_t,Nt::max()+1> counts_;

public:
    Counts() {
        for (size_t& c : counts_)
            c=-1;
        for (Nt nt : Nt::all)
            counts_[size_t(nt)] = 0;
    }
    Counts(const Counts<Nt4>& nt4counts);
    Counts(const Counts<Nt2>& nt2counts);

    void clear() {for (Nt nt : Nt::all) counts_[size_t(nt)]=0;}
    void increment(Nt nt) {++counts_[size_t(nt)];}
    void increment(Nt nt, size_t cnt) {counts_[size_t(nt)] += cnt;}

    // Get the count for a given nucleotide.
    size_t operator[] (Nt nt) const {return counts_[size_t(nt)];}
    const array<size_t,Nt::max()+1>& arr() const {return counts_;}

    size_t sum() const {return (*this)[Nt::a] + (*this)[Nt::c] + (*this)[Nt::g] + (*this)[Nt::t];}
    array<pair<size_t,Nt>,4> sorted() const;

    Counts& operator+= (const Counts& other)
        {for (Nt nt : Nt::all) counts_[size_t(nt)] += other.counts_[size_t(nt)]; return *this;}

    // Print the counts.
    template<typename Nt_> friend ostream& operator<< (ostream& os, const Counts<Nt_>& cnts);
};

//
// Inline definitions
// ==========
//

inline // (Note: Trivial, but this can't be defined in-class.)
Nt4::Nt4(const Nt2 nt2) : nt_(from_nt2[size_t(nt2)]) {}

template<> inline
Counts<Nt2>::Counts(const Counts<Nt2>& nt2counts) : counts_(nt2counts.counts_) {}
template<> inline
Counts<Nt2>::Counts(const Counts<Nt4>& nt4counts) : Counts() {
    for (Nt2 nt2 : Nt2::all) // We thus ignore Nt4::n.
        counts_[size_t(nt2)] = nt4counts[Nt4(nt2)];
}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt4>& nt4counts) : counts_(nt4counts.counts_) {}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt2>& nt2counts) : Counts() {
    for (Nt2 nt2 : Nt2::all)
        counts_[size_t(Nt4(nt2))] = nt2counts[nt2];
}

template<typename Nt>
array<pair<size_t,Nt>,4> Counts<Nt>::sorted() const {
    array<pair<size_t,Nt>,4> arr {{
        {(*this)[Nt::t], Nt::t},
        {(*this)[Nt::g], Nt::g},
        {(*this)[Nt::c], Nt::c},
        {(*this)[Nt::a], Nt::a}
    }};
    // Sort by decreasing value. Primarily on the count, secondarily on
    // the Nt4 value.
    std::sort(arr.rbegin(), arr.rend());
    return arr;
}

template<typename Nt>
ostream& operator<< (ostream& os, const Counts<Nt>& cnts) {
    bool first=true;
    for (Nt nt : Nt::all) {
        if (first)
            first = false;
        else
            os << " ";
        os << nt << ":" << cnts[nt];
    }
    return os;
}

#endif
