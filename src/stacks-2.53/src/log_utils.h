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

#ifndef __LOG_UTILS_H__
#define __LOG_UTILS_H__

#include <string>
#include <iostream>
#include <fstream>

#include "constants.h"
#include "utils.h"

int init_log(ostream &fh, int argc, char **argv);

// Returns e.g. "23.2%".
string as_percentage(double d);
inline string as_percentage(size_t n, size_t tot) {return as_percentage((double) n / tot);}

string to_string(const FileT& ft);

// ProgressMeter
// ==========
// Reports completeness to the given stream.
// In percentage mode reports 1-2-5-10-20-50-100%; in open-scale mode reports
// every 1eX, 2eX, 5eX starting at the passed parameter.
class ProgressMeter {
    ostream& os_;
    bool pct_;
    size_t n_max_;
    size_t n_done_;
    size_t next_;
    Timer timer_;

public:
    ProgressMeter(ostream& os, bool pct, size_t n_operations);
    ProgressMeter& operator+=(size_t n);
    ProgressMeter& operator++() {return operator+=(1);}
    void done();
};

// TeeBuf
// ==========
// A streambuf which tees to two streambufs.
// From http://wordaligned.org/articles/cpp-streambufs */
// This tee buffer has no actual buffer, so every character "overflows"
// directly into the teed buffers.
class TeeBuf: public streambuf {
public:

    TeeBuf(streambuf* sb1, streambuf* sb2)
        : sb1(sb1) , sb2(sb2)
        {}

private:
    streambuf* sb1;
    streambuf* sb2;

    virtual int overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            int r1 = sb1->sputc(c);
            int r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    // Sync both teed buffers.
    virtual int sync() {
        int r1 = sb1->pubsync();
        int r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
};

// LogAlterator
// ==========
class LogAlterator {
public:
    string log_path;
    string distribs_path;
    ofstream l; // The actual log file
    ostream  o; // Just stdout
    ostream  e; // Just stderr
    ofstream x; // The '.dist' file.

    // Construct the log alterator.
    // cout and cerr will also write to the log file (if quiet
    // is true, output to stdout and stderr is suppressed i.e. they only
    // write to the log file).
    LogAlterator(const string& prefix_path, bool distribs, bool quiet, int argc, char ** argv);

    // Upon destruction, restore cout and cerr.
    ~LogAlterator() {
        cout.rdbuf(o.rdbuf());
        cerr.rdbuf(e.rdbuf());
    }

private:
    TeeBuf lo_buf;
    TeeBuf le_buf;
};

#endif // __LOG_UTILS_H__
