// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cmath>
#include <ctime>
#include <iostream>
#include <utility>
#include <string>

#include <unistd.h>
#include <dirent.h>
#include <zlib.h>

#include "htslib/bgzf.h"
#include "constants.h"
#include "nucleotides.h"

char   reverse(char);
char  *rev_comp(const char *);
void   rev_comp_inplace(char*);
string rev_comp(const string&);
void   reverse_string(char *);
int    is_integer(const char *);
double is_double(const char *);
int    tokenize_string(const char *, vector<string> &);

double factorial(double);
double reduced_factorial(double, double);

double log_factorial(double);
double reduced_log_factorial(double, double);

inline
bool almost_equal(double x, double y) {
    const double precision = 1.0e-9;
    if (x == 0.0 && y == 0.0)
        return true;
    if (!std::isnormal(x) || !std::isnormal(y)) {
        stringstream ss;
        ss << "almost_equal: x=" << x << ", y=" << y;
        throw std::domain_error(ss.str());
    }
    return std::abs(x-y) <= precision * std::abs(std::min(x,y));
}

struct OnlineMeanVar {
    // Computes the mean and variance in a numerically stable way.
    // Uses the algorithm described in:
    // B. P. Welford. (1962) Note on a Method for Calculating Corrected Sums of
    // Squares and Products. Technometrics: 4(3), pp. 419-420.
    // Chan TF, Golub GH, LeVeque RJ (1979), "Updating Formulae and a Pairwise
    // Algorithm for Computing Sample Variances.", Technical Report STAN-CS-79-773,
    // Dpt. Comp. Sci., Stanford University.
    double n_;
    double mean_;
    double M2_;
public:
    OnlineMeanVar() : n_(0.0), mean_(0.0), M2_(0.0) {}

    size_t n()     const {return std::round(n_);}
    double n_dbl() const {return n_;}
    double mean()  const {return mean_;}
    double var_p() const {return n_ > 0 ? M2_/n_ : NAN;}
    double var_s() const {return n_ > 1 ? M2_/(n_-1) : NAN;}
    double sd_p()  const {return std::sqrt(var_p());}
    double sd_s()  const {return std::sqrt(var_s());}

    void increment(double x) {
        n_ += 1.0;
        double delta1 = x - mean_;
        mean_ += delta1 / n_;
        double delta2 = x - mean_;
        M2_ += delta1 * delta2;
    }

    OnlineMeanVar& operator+=(const OnlineMeanVar& other) {
        double delta = other.mean_ - mean_;
        double n_tot = n_ + other.n_;
        double w_other = other.n_ / n_tot;
        mean_ += delta * w_other;
        M2_ += other.M2_ + delta * delta * n_ * w_other;
        n_ = n_tot;
        return *this;
    }
};

//
// GtLiks: A class to store the likelihoods of SNP genotypes.
//
class GtLiks {
    array<double,10> lnliks_; // {AA,AC,CC,AG,CG,GG,AT,CT,GT,TT} similar to VCF.
public:
    GtLiks() : lnliks_{{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}} {}
    double at(Nt2 n1, Nt2 n2) const {return at(gt_index(n1,n2));}
    bool has_lik(Nt2 n1, Nt2 n2) const {return has_lik(gt_index(n1,n2));}
    void set(Nt2 n1, Nt2 n2, double lnl) {set(gt_index(n1, n2), lnl);}

    double at(size_t gt) const {assert(has_lik(gt)); return lnliks_[gt];}
    bool has_lik(size_t gt) const {return lnliks_[gt] != 1.0;}
    void set(size_t gt, double lnl) {assert(!std::isnan(lnl)); assert(lnl<=0.0); assert(!has_lik(gt)); lnliks_[gt] = lnl;}

    static size_t gt_index(Nt2 n1, Nt2 n2) {
        if(n1<n2)
            return size_t(n1) + (size_t(n2)*(size_t(n2)+1)) / 2;
        else
            return size_t(n2) + (size_t(n1)*(size_t(n1)+1)) / 2;
    }

    // For debugging.
    friend ostream& operator<<(ostream& os, const GtLiks& liks);
};

//
// Comparison functions for the STL sort routine
//
bool compare_ints(int, int);
bool compare_pair(pair<char, int>, pair<char, int>);
bool compare_pair_intint(pair<int, int>, pair<int, int>);
bool compare_pair_intdouble(pair<int, double>, pair<int, double>);
bool compare_pair_stringint(pair<string, int>, pair<string, int>);
bool compare_pair_haplotype(const pair<string, double>&, const pair<string, double>&);
bool compare_pair_haplotype_rev(const pair<string, double>&, const pair<string, double>&);
bool compare_str_len(const string&, const string&);
struct LessCStrs
{
   bool operator() (const char* str1, const char* str2) const {return strcmp(str1, str2) < 0 ;}
} ;
struct int_increasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs < rhs;
    }
};

struct int_decreasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs > rhs;
    }
};

//
// Remove elements from a vector.
//
template<typename T, typename Predicate>
void stacks_erase_if(vector<T>& v, Predicate p)
{
    v.erase(std::remove_if(v.begin(), v.end(), p), v.end());
}

//
// Join a range of elements into a stream.
//
template<typename IterableT, typename SepT>
void join(IterableT elements, const SepT& sep, ostream& os) {
    auto first = elements.begin();
    if (first != elements.end()) {
        os << *first;
        ++first;
        while (first != elements.end()) {
            os << sep << *first;
            ++first;
        }
    }
}

//
// Routines to check that files are open.
//
inline
void check_open (std::ifstream& fs, const string& path) {
    if (!fs.is_open()) {
        cerr << "Error: Failed to open '" << path << "' for reading.\n";
        throw exception();
    }
    fs.exceptions(fs.exceptions() | ios::badbit);
}
inline
void check_open (std::ofstream& fs, const string& path) {
    if (!fs.is_open()) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }
    fs.exceptions(fs.exceptions() | ios::badbit);
}
inline
void check_open (const gzFile fs, const string& path)
    {if (fs == NULL) {cerr << "Error: Failed to gz-open file '" << path << "'.\n"; throw exception();}}
inline
void check_open (const BGZF* fs, const string& path)
    {if (fs == NULL) {cerr << "Error: Failed to bgzf-open file '" << path << "'.\n"; throw exception();}}

//
// Check that a directory exists or try to create it.
//
void check_or_mk_dir(const string& path);

//
// Chronometer.
//
class Timer {
    double elapsed_;
    double consumed_;
    double start_;

public:
    Timer() : elapsed_(0.0), consumed_(0.0), start_(0.0) {restart();}
    void   restart() {start_=gettm(); consumed_+=2.0*(gettm()-start_);}
    void   update()  {double now=gettm(); elapsed_+=now-start_; start_=now; consumed_+=2.0*(gettm()-now);}
    double elapsed() const {return elapsed_;}
    double consumed() const {return consumed_;}

    Timer& operator+=(const Timer& other) {elapsed_+=other.elapsed_; consumed_+=other.consumed_; return *this;}

private:
    double gettm() {
        #if HAVE_CLOCK_GETTIME && defined _POSIX_MONOTONIC_CLOCK && _POSIX_MONOTONIC_CLOCK >= 0
        struct timespec ts;
        if(clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
            return ts.tv_sec + ts.tv_nsec / 1.0e9;
        else
            return 0.0;
        #else
        return 0.0;
        #endif
    }
};

//
// Class to read lines from a plain text or compressed file indifferently.
//
class VersatileLineReader {
    string path_;
    size_t line_number_;
    bool is_gzipped_;

    ifstream ifs_;
    string ifsbuffer_;

    gzFile gzfile_;
    char* gzbuffer_;
    size_t gzbuffer_size_;
    size_t gzline_len_;
    static const size_t gzbuffer_init_size = 65536;

public:
    VersatileLineReader(const string& path);
    VersatileLineReader();
    ~VersatileLineReader();

    int open(string &);

    //
    // Reads one line from the file, removing the trailing '\n' (and '\r', if any).
    // Returns false on EOF, or throws an exception if the file doesn't end with a newline.
    // e.g.:
    // const char* line; size_t len; while (file.getline(line, len)) {...}
    //
    bool getline(const char*& line, size_t& len);

    const string& path() const {return path_;}
    size_t line_number() const {return line_number_;} // 1-based.
};

class VersatileWriter {
    const string path_;
    bool is_gzipped_;
    ofstream ofs_;
    gzFile gzfile_;

    void gzputs_(const char* s);

public:
    VersatileWriter(const string& path);
    ~VersatileWriter() {gzclose(gzfile_);}

    const string& path() const {return path_;}
    void close();

    friend VersatileWriter& operator<< (VersatileWriter& w, char c);
    friend VersatileWriter& operator<< (VersatileWriter& w, const char* s);
    friend VersatileWriter& operator<< (VersatileWriter& w, const string& s);
    friend VersatileWriter& operator<< (VersatileWriter& w, int i);
    friend VersatileWriter& operator<< (VersatileWriter& w, long i);
    friend VersatileWriter& operator<< (VersatileWriter& w, size_t i);
};

//
// Wrapper for directory parsing functions.
// e.g. for(DirIterator e (path); e; ++e) {...}
//
class DirIterator {
    DIR* dir;
    struct dirent* entry;
public:
    DirIterator(const string& dir_path) : dir(NULL), entry(NULL) {
        dir = opendir(dir_path.c_str());
        if (dir == NULL) {
            cerr << "Error: Unable to open directory '" << dir_path << "' for reading.\n";
            throw exception();
        }
        entry = readdir(dir);
    }
    ~DirIterator() {closedir(dir);}

    const char* name() const {return entry->d_name;}

    operator bool() const {return entry!=NULL;}
    DirIterator& operator++() {entry = readdir(dir); return *this;}
    dirent* operator*() {return entry;}
};

// strip_read_number
// ----------
// Given a read name, removes the trailing /1, /2, _1 or _1.
inline
void strip_read_number(string& read_name) {
    if (read_name.size() >= 2) {
        // Check for ".../1" or ".../2"
        if ((*read_name.rbegin() == '1' || *read_name.rbegin() == '2')
            && (*++read_name.rbegin() == '/' || *++read_name.rbegin() == '_')
        ) {
            // Remove the suffix & return.
            read_name.resize(read_name.size()-2);
            return;
        }
        // Check for "... 1:..." or "... 2:..."
        const char *p, *q;
        if ((p = q = strchr(read_name.c_str(), ' ')) != NULL) {
            if ((*++q == '1' || *q == '2') && *++q == ':') {
                read_name.resize(p - read_name.c_str());
                return;
            }
        }
    }
    // Unexpected suffix.
    cerr << "Error: Unrecognized paired-end read name format, at '" << read_name << "'.\n";
    throw exception();
}

inline
void VersatileWriter::close() {
    if (is_gzipped_) {
        if (gzclose(gzfile_) != Z_OK)
            throw ios::failure("gzclose");
        gzfile_ = NULL;
    } else {
        ofs_.close();
    }
}

inline
void gzclose_throwing(gzFile f) {
    if (gzclose(f) != Z_OK)
        throw ios::failure("gzclose");
}

inline
VersatileWriter& operator<< (VersatileWriter& w, char c) {
    if (w.is_gzipped_) {
        if (gzputc(w.gzfile_, c) == -1)
            throw ios::failure("gzputc");
    } else {
        w.ofs_ << c;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, const char* s) {
    if (w.is_gzipped_)
        w.gzputs_(s);
    else
        w.ofs_ << s;
    return w;
}

inline
void VersatileWriter::gzputs_(const char* s) {
    if (gzputs(gzfile_, s) == -1)
        throw ios::failure("gzputs");
}

inline
void gzputs_throwing(gzFile f, const char* s) {
    if (gzputs(f, s) == -1)
        throw ios::failure("gzputs");
}

inline
VersatileWriter& operator<< (VersatileWriter& w, const string& s) {
    if (w.is_gzipped_) {
        if (!s.empty() && gzwrite(w.gzfile_, s.c_str(), s.length()) <= 0)
            throw ios::failure("gzwrite");
    } else {
        w.ofs_ << s;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, int i) {
    if (w.is_gzipped_) {
        char buf[16];
        sprintf(buf, "%d", i);
        w.gzputs_(buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, long i) {
    if (w.is_gzipped_) {
        char buf[32];
        sprintf(buf, "%ld", i);
        w.gzputs_(buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, size_t i) {
    if (w.is_gzipped_) {
        char buf[32];
        sprintf(buf, "%zu", i);
        w.gzputs_(buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

#endif // __UTILS_H__
