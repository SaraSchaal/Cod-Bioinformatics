// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2015, Julian Catchen <jcatchen@illinois.edu>
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

#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "constants.h"
#include "utils.h"
#include "log_utils.h"

using namespace std;

int
init_log(ostream &fh, int argc, char **argv)
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    //
    // Write the command line that was executed.
    //
    string exe (argv[0]);
    fh << exe.substr(exe.find_last_of('/')+1) << " v" << VERSION << ", executed " << date;
    #ifdef HAVE_LIBZ
    fh << " (zlib-" << zlibVersion() << ")";
    #endif
    fh << "\n";
    for (int i = 0; i < argc; i++) {
        fh << argv[i];
        if (i < argc - 1) fh << " ";
    }
    fh << "\n";

    fh << flush;
    return 0;
}

class Pct {
    double d_;
public:

    friend ostream& operator<< (ostream& os, const Pct& p) {
        ostream os2 (os.rdbuf());
        os2 << std::fixed << setprecision(1);
        os2 << p.d_ * 100 << "%";
        return os;
    }
};

string as_percentage(double d) {
    stringstream ss;
    ss << std::fixed << setprecision(1);
    ss << d*100 << "%";
    return ss.str();
}

string to_string(const FileT& ft) {
    if (ft == FileT::unknown)
        return "unknown";
    if (ft == FileT::sql)
        return "sql";
    if (ft == FileT::gzsql)
        return "gzsql";
    if (ft == FileT::fasta)
        return "fasta";
    if (ft == FileT::gzfasta)
        return "gzfasta";
    if (ft == FileT::fastq)
        return "fastq";
    if (ft == FileT::gzfastq)
        return "gzfastq";
    if (ft == FileT::bowtie)
        return "bowtie";
    if (ft == FileT::sam)
        return "sam";
    if (ft == FileT::bam)
        return "bam";
    if (ft == FileT::tsv)
        return "tsv";
    if (ft == FileT::bustard)
        return "bustard";
    if (ft == FileT::phase)
        return "phase";
    if (ft == FileT::fastphase)
        return "fastphase";
    if (ft == FileT::beagle)
        return "beagle";
    DOES_NOT_HAPPEN;
    return "?!";
}

ProgressMeter::ProgressMeter(
        ostream& os,
        bool pct,
        size_t n_operations
): os_(os),
   pct_(pct),
   n_max_(pct_ ? n_operations : SIZE_MAX),
   n_done_(0),
   next_(pct_ ? n_max_*0.01 : n_operations)
{}

ProgressMeter& ProgressMeter::operator+=(size_t n)
{
    assert(n_done_ < n_max_);
    n_done_ += n;

    if (n_done_ >= next_) {
        #ifdef DEBUG
        timer_.update();
        os_ << '(' << (size_t) timer_.elapsed() << "s) ";
        #endif
        if (pct_) {
            if (n_done_ >= size_t(n_max_ * 0.5)) {
                os_ << "50%...\n";
                next_ = SIZE_MAX;
            } else if (n_done_ >= size_t(n_max_ * 0.2)) {
                os_ << "20%...\n";
                next_ = n_max_ * 0.5;
            } else if (n_done_ >= size_t(n_max_ * 0.1)) {
                os_ << "10%...\n";
                next_ = n_max_ * 0.2;
            } else if (n_done_ >= size_t(n_max_ * 0.05)) {
                os_ << "5%...\n";
                next_ = n_max_ * 0.1;
            } else if (n_done_ >= size_t(n_max_ * 0.02)) {
                os_ << "2%...\n";
                next_ = n_max_ * 0.05;
            } else {
                os_ << "1%...\n";
                next_ = n_max_ * 0.02;
            }
        } else {
            if (n_done_ >= 1000)
                os_ << n_done_/1000 << 'K';
            else
                os_ << n_done_;
            os_ << "...\n";
            size_t scale = 1;
            while (scale <= next_)
                scale *= 10;
            if (next_ < scale / 5)
                next_ = scale / 5;
            else if (next_ < scale / 2)
                next_ = scale / 2;
            else
                next_ = scale;
        }
        os_ << flush;
    }

    return *this;
}

void ProgressMeter::done()
{
    #ifdef DEBUG
    timer_.update();
    os_ << "(" << (size_t) timer_.elapsed() << "s) ";
    #endif
    if (pct_) {
        assert(n_done_ == n_max_);
        os_ << "100%\n";
    } else {
        os_ << "done.\n";
    }
}

LogAlterator::LogAlterator(const string& prefix, bool distribs, bool quiet, int argc, char** argv)
    : log_path(prefix + ".log")
    , distribs_path(prefix + ".log.distribs")
    , l(log_path)
    , o(cout.rdbuf())
    , e(cerr.rdbuf())
    , x()
    , lo_buf(cout.rdbuf(), l.rdbuf())
    , le_buf(cerr.rdbuf(), l.rdbuf())
{
    check_open(l, log_path);
    init_log(l, argc, argv);
    if (distribs) {
        x.open(distribs_path);
        check_open(x, distribs_path);
        x << "# Note: Individual distributions can be extracted using the `stacks-dist-extract` utility.\n"
          << "#       e.g. `stacks-dist-extract "
          << distribs_path.c_str() + (distribs_path.find_last_of('/') + 1)
          <<" dist_name`\n";
    }

    if (quiet) {
        // Use the fstream buffer only, and not at all stdout and stderr.
        cout.rdbuf(l.rdbuf());
    } else {
        cout << "Logging to '" << log_path << "'.\n";
        cout.rdbuf(&lo_buf);
        if (distribs)
            cout << "Locus/sample distributions will be written to '" << distribs_path << "'.\n";
    }
    cerr.rdbuf(&le_buf);
}
