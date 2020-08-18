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

//
// write.cc -- common routines for writing FASTA/FASTQ records to a file..
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "write.h"

int
write_fasta(ofstream *fh, RawRead *href, bool overhang) {
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
        *fh <<
            ">" << href->run <<
            "_" << href->lane <<
            "_" << tile <<
            "_" << href->x <<
            "_" << href->y <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n";
    else
        *fh <<
            ">" << href->machine <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fasta(gzFile *fh, RawRead *href, bool overhang) {
    stringstream sstr;
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
        sstr <<
            ">" << href->run <<
            "_" << href->lane <<
            "_" << tile <<
            "_" << href->x <<
            "_" << href->y <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n";
    else
        sstr <<
            ">" << href->machine <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fasta(ofstream *fh, Seq *href) {
    *fh <<
        ">" <<
        href->id  << "\n" <<
        href->seq << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fasta(gzFile *fh, Seq *href) {
    stringstream sstr;

    sstr <<
        ">" <<
        href->id  << "\n" <<
        href->seq << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fastq(ofstream *fh, RawRead *href, bool overhang) {
    //
    // Write the sequence and quality scores in FASTQ format.
    //
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
        *fh <<
            "@" << href->run <<
            "_" << href->lane <<
            "_" << tile <<
            "_" << href->x <<
            "_" << href->y <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n" <<
            "+\n" <<
            href->phred + offset << "\n";
    else
        *fh <<
            "@" << href->machine <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n" <<
            "+\n" <<
            href->phred + offset << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fastq(gzFile *fh, RawRead *href, bool overhang) {
    //
    // Write the sequence and quality scores in FASTQ format.
    //
    stringstream sstr;
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
        sstr <<
            "@" << href->run <<
            "_" << href->lane <<
            "_" << tile <<
            "_" << href->x <<
            "_" << href->y <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n" <<
            "+\n" <<
            href->phred + offset << "\n";
    else
        sstr <<
            "@" << href->machine <<
            "/" << href->read << "\n" <<
            href->seq + offset << "\n" <<
            "+\n" <<
            href->phred + offset << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fastq(ofstream *fh, Seq *href, int offset) {
    *fh <<
        "@" << href->id     << "\n" <<
        href->seq + offset  << "\n" <<
        "+\n" <<
        href->qual + offset << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fastq(gzFile *fh, Seq *href, int offset) {
    stringstream sstr;
    sstr <<
        "@" << href->id     << "\n" <<
        href->seq + offset  << "\n" <<
        "+\n" <<
        href->qual + offset << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fasta(ofstream *fh, Seq *href, int offset) {
    *fh <<
        ">" <<
        href->id << "\n" <<
        href->seq + offset << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fasta(gzFile *fh, Seq *href, int offset) {
    stringstream sstr;
    sstr <<
        ">" <<
        href->id << "\n" <<
        href->seq + offset << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fastq(ofstream *fh, Seq *href) {
    *fh <<
        "@" << href->id << "\n" <<
        href->seq << "\n" <<
        "+\n" <<
        href->qual << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fastq(gzFile *fh, Seq *href) {
    stringstream sstr;
    sstr <<
        "@" << href->id << "\n" <<
        href->seq << "\n" <<
        "+\n" <<
        href->qual << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fastq(ofstream *fh, Seq *href, string msg) {
    *fh <<
        "@" << href->id << "|" << msg << "\n" <<
        href->seq << "\n" <<
        "+\n" <<
        href->qual << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fastq(gzFile *fh, Seq *href, string msg) {
    stringstream sstr;
    sstr <<
        "@" << href->id << "|" << msg << "\n" <<
        href->seq << "\n" <<
        "+\n" <<
        href->qual << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fasta(ofstream *fh, Seq *href, string msg) {
    *fh <<
        ">" <<
        href->id  << "|" << msg << "\n" <<
        href->seq << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fasta(gzFile *fh, Seq *href, string msg) {
    stringstream sstr;
    sstr <<
        ">" <<
        href->id  << "|" << msg << "\n" <<
        href->seq << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fasta(ofstream *fh, Seq *href, RawRead *r) {
    *fh        << ">"
        << href->id << "\n"
        << r->seq + r->inline_bc_len << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fasta(gzFile *fh, Seq *href, RawRead *r) {
    stringstream sstr;
    sstr << ">"
         << href->id << "\n"
         << r->seq + r->inline_bc_len << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}

int
write_fastq(ofstream *fh, Seq *href, RawRead *r) {
    *fh << "@" << href->id << "\n"
        << r->seq   + r->inline_bc_len << "\n"
        << "+\n"
        << r->phred + r->inline_bc_len << "\n";

    if (fh->fail()) return -1;

    return 1;
}

int
write_fastq(gzFile *fh, Seq *href, RawRead *r) {
    stringstream sstr;
    sstr << "@" << href->id << "\n"
         << r->seq   + r->inline_bc_len << "\n"
         << "+\n"
         << r->phred + r->inline_bc_len << "\n";

    int res = gzputs(*fh, sstr.str().c_str());

    return res;
}
