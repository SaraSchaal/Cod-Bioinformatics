#include "DNASeq4.h"

DNASeq4::DNASeq4(const char* s, size_t len) : l_(len), v_() {
    reserve(l_);
    for (size_t i=0; i<len; i+=2)
        v_.push_back(DiNuc(s[i], s[i+1])); //n.b. `s` is null-terminated
}

string DNASeq4::str() const {
    string s;
    s.reserve(l_);
    for (auto nt=begin(); nt!=end(); ++nt)
        s.push_back(char(*nt));
    return s;
}

DNASeq4 DNASeq4::rev_compl() const {
    DNASeq4 rev;
    rev.l_ = l_;
    rev.reserve(rev.l_);

    for (iterator nt=end(); nt!=begin();) {
        --nt;
        Nt4 first = (*nt).rev_compl();
        if (nt != begin()) {
            --nt;
            Nt4 second = (*nt).rev_compl();
            rev.v_.push_back(DiNuc(first, second));
        } else {
            rev.v_.push_back(DiNuc((*nt).rev_compl(), Nt4::$));
            break;
        }
    }

    return rev;
}

void DNASeq4::resize(size_t len) {
    if (l_ < len) {
        while (l_ != len) {
            push_back(Nt4::n);
            ++l_;
        }
    } else if (l_ > len) {
        l_ = len;
        v_.resize(l_/2 + l_%2);
        if (l_%2)
            v_.back().second(Nt4::$);
    }
}

void DNASeq4::append(iterator first, iterator past) {

    if (!(first != past))
        return;

    if (l_%2==1) {
        v_.back().second(*first);
        ++l_;
        ++first;
    }

    l_ += past - first;
    reserve(l_);

    while(first != past) {
        Nt4 prev = *first;
        ++first;
        if (first != past) {
            v_.push_back(DiNuc(prev, *first));
            ++first;
        } else {
            v_.push_back(DiNuc(prev, Nt4::$));
            break;
        }
    }
}

void DNASeq4::remove_Ns() {
    auto itr = begin();
    while (itr != end() && *itr != Nt4::n)
        ++itr;
    if (itr == end())
        // Sequence didn't contain N's.
        return;
    DNASeq4 seq;
    auto start = begin();
    while (start != end()) {
        while (itr != end() && *itr != Nt4::n)
            ++itr;
        seq.append(start, itr);
        while (itr != end() && *itr == Nt4::n)
            ++itr;
        start = itr;
    }
    *this = move(seq);
}

void DNASeq4::shift_Ns_towards_the_end() {
    vector<DiNuc>::iterator v_itr = v_.begin();
    iterator first = begin();
    while(first != end()) {
        while (first != end() && *first == Nt4::n)
            ++first;
        if (first == end())
            break;
        assert(*first != Nt4::n);

        iterator second = first;
        ++second;
        while (second != end() && *second == Nt4::n)
            ++second;
        if (second == end())
            break;
        assert(*second != Nt4::n);

        assert(v_itr != v_.end());
        *v_itr = DiNuc(*first, *second);
        ++v_itr;

        first = second;
        ++first;
    }
    if (first != end()) {
        // Loop ended while scanning for a `second`, and nucleotide `*first`
        // hasn't been written yet.
        assert(v_itr != v_.end());
        *v_itr = DiNuc(*first, Nt4::n);
        ++v_itr;
    }
    while (v_itr != v_.end()) {
        *v_itr = DiNuc(Nt4::n, Nt4::n);
        ++v_itr;
    }
    if (l_%2)
        v_.back().second(Nt4::$);
}
