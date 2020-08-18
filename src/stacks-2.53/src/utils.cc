// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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
// utils.cc -- common routines needed in multiple object files.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include <regex>

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "utils.h"

using namespace std;

char reverse(char c) {
    switch (c) {
    case 'A':
    case 'a':
        return 'T';
        break;
    case 'C':
    case 'c':
        return 'G';
        break;
    case 'G':
    case 'g':
        return 'C';
        break;
    case 'T':
    case 't':
        return 'A';
        break;
    case 'N':
    case 'n':
    case '.':
        return 'N';
        break;
    case '-':
    default:
        return '-';
        break;
    }

    return 'N';
}

string
rev_comp(const string& seq)
{
    string s;
    s.reserve(seq.length());
    for (char c: seq) {
        switch (c) {
        case 'A':
        case 'a':
            s.push_back('T');
            break;
        case 'C':
        case 'c':
            s.push_back('G');
            break;
        case 'G':
        case 'g':
            s.push_back('C');
            break;
        case 'T':
        case 't':
            s.push_back('A');
            break;
        default:
            s.push_back('N');
            break;
        }
    }

    return s;
}

char *
rev_comp(const char *seq)
{
    int len   = strlen(seq);
    int j     = 0;
    char *com = new char[len + 1]; // *** new
    const char *p;

    for (p = seq + len - 1; p >= seq; p--) {
        switch (*p) {
        case 'A':
        case 'a':
            com[j] = 'T';
            break;
        case 'C':
        case 'c':
            com[j] = 'G';
            break;
        case 'G':
        case 'g':
            com[j] = 'C';
            break;
        case 'T':
        case 't':
            com[j] = 'A';
            break;
        case 'N':
        case 'n':
        case '.':
            com[j] = 'N';
            break;
        }
        j++;
    }
    com[len] = '\0';

    return com;
}

void rev_comp_inplace(char* seq) {
    assert(seq);
    char* p = seq;
    while(*p) {
        switch (*p) {
        case 'A':
        case 'a':
            *p = 'T';
            break;
        case 'C':
        case 'c':
            *p = 'G';
            break;
        case 'G':
        case 'g':
            *p = 'C';
            break;
        case 'T':
        case 't':
            *p = 'A';
            break;
        case 'N':
        case 'n':
        case '.':
            *p = 'N';
            break;
        }
        ++p;
    }
    std::reverse(seq, p);
}

void
reverse_string(char *seq)
{
    int len = strlen(seq);
    char *p = seq;
    char *q = seq + len - 1;
    char  tmp;

    while (q > p) {
        tmp = *q;
        *q  = *p;
        *p  = tmp;
        q--;
        p++;
    }

    return;
}

int
is_integer(const char *str)
{
    //
    // Adapted from the strtol manpage.
    //
    char *endptr;

    // To distinguish success/failure after call
    errno = 0;
    long val = strtol(str, &endptr, 10);

    //
    // Check for various possible errors
    //
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0)) {
        return -1;
    }

    if (endptr == str || *endptr != '\0')
        return -1;

    return (int) val;
}

double
is_double(const char *str)
{
    //
    // Adapted from the strtol manpage.
    //
    char *endptr;

    // To distinguish success/failure after call
    errno = 0;
    double val = strtod(str, &endptr);

    //
    // Check for various possible errors
    //
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0)) {
        return -1;
    }

    if (endptr == str || *endptr != '\0')
        return -1;

    return val;
}

double
factorial(double n)
{
    double fact = 1;

    for (double i = n; i > 1; i--)
        fact *= i;

    return fact;
}

double
reduced_factorial(double n, double d)
{
    double f = n - d;

    if (f < 0)
        return 0;
    else if (f == 0)
        return 1;
    else if (f == 1)
        return n;

    f = n;
    n--;
    while (n > d) {
        f *= n;
        n--;
    }

    return f;
}

double
log_factorial(double n)
{
    double fact = 0;

    for (double i = n; i > 1; i--)
        fact += log(i);

    return fact;
}

double
reduced_log_factorial(double n, double d)
{
    double f = n - d;

    if (f < 0)
        return 0;
    else if (f == 0)
        return 0;
    else if (f == 1)
        return log(n);

    f = log(n);
    n--;
    while (n > d) {
        f += log(n);
        n--;
    }

    return f;
}

ostream& operator<<(ostream& os, const GtLiks& liks) {
    bool first = true;
    for (auto nt2=Nt2::all.begin(); nt2!=Nt2::all.end(); ++nt2) {
        for (auto nt1=Nt2::all.begin(); ; ++nt1) {
            if (first)
                first = false;
            else
                os << ", ";
            if (liks.has_lik(*nt1, *nt2))
                cout << *nt1 << *nt2 << ":" << liks.at(*nt2, *nt1);
            else
                cout << char(std::tolower(char(*nt1))) << char(std::tolower(char(*nt2)));
            if (nt1==nt2)
                break;
        }
    }
    return os;
}

//
// Tokenize a string.
//
int
tokenize_string(const char *buf, vector<string> &tokens)
{
    const char *p, *q;

    p = buf;

    do {
        q = p;

        while (*q != '\0' && *q != ' ' && *q != '\t') q++;
        tokens.push_back(string(p, q - p + 1));

        p = q + 1;

    } while (*q != '\0');

    return 0;
}


bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}

bool compare_pair_intint(pair<int, int> a, pair<int, int> b) {
    return (a.second > b.second);
}

bool compare_pair_intdouble(pair<int, double> a, pair<int, double> b) {
    return (a.second < b.second);
}

bool compare_pair_stringint(pair<string, int> a, pair<string, int> b) {
    return (a.second < b.second);
}

bool compare_ints(int a, int b) {
    return (a > b);
}

bool compare_pair_haplotype(const pair<string, double>& a, const pair<string, double>& b) {
    return (a.second > b.second);
}

bool compare_pair_haplotype_rev(const pair<string, double>& a, const pair<string, double>& b) {
    return (a.second < b.second);
}

bool compare_str_len(const string& a, const string& b) {
    return (a.length() < b.length());
}

void check_or_mk_dir(const string& path) {
    string path_stripped;
    const string* path_p = &path;
    if (path.back() == '/') {
        path_stripped = path;
        path_stripped.pop_back();
        path_p = &path_stripped;
    }

    struct stat s;
    if (stat(path_p->c_str(), &s) == 0) {
        //
        // Path exists, check that it is a directory
        //
        if (!S_ISDIR(s.st_mode)) {
            cerr << "Error: '" << *path_p << "' is not a directory.\n";
            throw exception();
        }
    } else if (mkdir(path_p->c_str(), ACCESSPERMS) != 0) {
        //
        // Failed to create the directory.
        //
        cerr << "Error: Failed to create directory '" << *path_p << "'.\n";
        throw exception();
    }
}

VersatileLineReader::VersatileLineReader()
    : path_(), line_number_(0), is_gzipped_(false),
      ifs_(), ifsbuffer_(),
      gzfile_(NULL), gzbuffer_(NULL), gzbuffer_size_(0), gzline_len_(0) {}

VersatileLineReader::VersatileLineReader(const string& path)
    : path_(path), line_number_(0), is_gzipped_(false),
      ifs_(), ifsbuffer_(),
      gzfile_(NULL), gzbuffer_(NULL), gzbuffer_size_(0), gzline_len_(0)
{
    std::smatch m;
    std::regex_search(path_, m, std::regex("\\.[Gg][Zz]$"));
    is_gzipped_ = !m.empty();

    if (!is_gzipped_) {
        ifs_.open(path_);
        check_open(ifs_, path_);
    } else {
        gzfile_ = gzopen(path_.c_str(), "rb");
        check_open(gzfile_, path_);
        gzbuffer_size_ = gzbuffer_init_size;
        gzbuffer_ = new char[gzbuffer_size_];
    }
}

VersatileLineReader::~VersatileLineReader() {
    if (is_gzipped_) {
        gzclose(gzfile_);
        delete[] gzbuffer_;
    }
}

int
VersatileLineReader::open(string &path)
{
    this->path_ = path;

    std::smatch m;
    std::regex_search(this->path_, m, std::regex("\\.[Gg][Zz]$"));
    this->is_gzipped_ = !m.empty();

    if (path.compare(path.length()-strlen("catalog.calls"), string::npos, "catalog.calls") == 0)
        is_gzipped_ = true;

    if (!this->is_gzipped_) {
        this->ifs_.open(this->path_);
        check_open(this->ifs_, this->path_);
    } else {
        this->gzfile_ = gzopen(this->path_.c_str(), "rb");
        check_open(this->gzfile_, this->path_);
        this->gzbuffer_size_ = gzbuffer_init_size;
        this->gzbuffer_ = new char[this->gzbuffer_size_];
    }

    return 0;
}

bool VersatileLineReader::getline(const char*& line, size_t& len) {
    auto truncated = [this]() {
        cerr << "Error: While reading '" << path_ << "' (file may be truncated).\n";
        throw exception();
    };

    if (!is_gzipped_) {
        if (!std::getline(ifs_, ifsbuffer_))
            return false;
        if (ifs_.eof())
            // Doesn't end with an '\n'.
            truncated();
        if (ifsbuffer_.back() == '\r')
            // Remove the '\r'.
            ifsbuffer_.pop_back();

        line = ifsbuffer_.c_str();
        len = ifsbuffer_.length();

    } else {
        if (!gzgets(gzfile_, gzbuffer_, gzbuffer_size_))
            return false;

        // Check the contents of the buffer.
        gzline_len_ = strlen(gzbuffer_);
        while (gzbuffer_[gzline_len_ - 1] != '\n') {
            if (gzline_len_ < gzbuffer_size_ - 1) {
                // The buffer was long enough but file didn't end with a newline.
                truncated();
            } else {
                // The buffer wasn't long enough.
                assert(gzline_len_ == gzbuffer_size_ - 1);

                // Get a new, wider buffer & copy what we've already read.
                char* old_buffer = gzbuffer_;
                gzbuffer_size_ *= 2;
                gzbuffer_ = new char[gzbuffer_size_];
                memcpy(gzbuffer_, old_buffer, gzline_len_ + 1);
                delete[] old_buffer;

                // Continue reading the line.
                char* start = gzbuffer_ + gzline_len_;
                assert(*start == '\0');
                if (!gzgets(gzfile_, start, gzbuffer_size_ - gzline_len_))
                    // EOF; this means that the last character read by the previous
                    // call was the last of the file. And it wasn't a newline.
                    truncated();
                gzline_len_ += strlen(start);
            }
        }
        // Remove the '\n'.
        gzbuffer_[gzline_len_ - 1] = '\0';
        --gzline_len_;
        if (gzbuffer_[gzline_len_ - 1] == '\r') {
            // Remove the '\r'.
            gzbuffer_[gzline_len_ - 1] = '\0';
            --gzline_len_;
        }

        line = gzbuffer_;
        len = gzline_len_;
    }

    ++line_number_;
    return true;
}

VersatileWriter::VersatileWriter(const string& path)
: path_(path),
  is_gzipped_(false),
  ofs_(),
  gzfile_()
{
    std::smatch m;
    std::regex_search(path, m, std::regex("\\.[Gg][Zz]$"));
    is_gzipped_ = !m.empty();

    if (path.compare(path.length()-strlen("catalog.calls"), string::npos, "catalog.calls") == 0)
        is_gzipped_ = true;

    if (!is_gzipped_) {
        ofs_.open(path_);
        check_open(ofs_, path_);
    } else {
        gzfile_ = gzopen(path_.c_str(), "wb");
        check_open(gzfile_, path_);
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gzfile_, libz_buffer_size);
        #endif
    }
}
