// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2015-2019, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __BAMUNALIGNEDI_H__
#define __BAMUNALIGNEDI_H__

//
// Code to parse binary BAM format. This format is created for
// reads that have NOT been aligned to a reference genome.
//

#include "input.h"
#include "htslib/sam.h"

class BamUnAln: public Input {
    htsFile   *bam_fh;
    bam_hdr_t *bamh;
    bam1_t    *aln;

    map<uint, string> chrs;

    int parse_header();

 public:
    BamUnAln(const char *path) : Input() {
        this->path   = string(path);
        this->bam_fh = hts_open(path, "r");
        this->aln    = bam_init1();

        this->parse_header();
    };
    BamUnAln(string path) : Input() {
        this->path   = path;
        this->bam_fh = hts_open(path.c_str(), "r");
        this->aln    = bam_init1();

        this->parse_header();
    };
    ~BamUnAln() {
        hts_close(this->bam_fh);
        bam_hdr_destroy(this->bamh);
        bam_destroy1(this->aln);
    };
    Seq *next_seq();
    int  next_seq(Seq &);
};

int
BamUnAln::parse_header()
{
    this->bamh = bam_hdr_init();
    this->bamh = sam_hdr_read(this->bam_fh);

    for (uint j = 0; j < (uint) this->bamh->n_targets; j++) {
        //
        // Record the mapping from integer ID to chromosome name that we will see in BAM records.
        //
        this->chrs[j] = string(this->bamh->target_name[j]);
    }

    return 0;
}

Seq *
BamUnAln::next_seq()
{
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int
BamUnAln::next_seq(Seq& s)
{
    int bytes_read = 0;

    //
    // Read a record from the file and place it in a Seq object.
    //
    bytes_read = sam_read1(this->bam_fh, this->bamh, this->aln);

    if (bytes_read <= 0)
        return 0;

    //
    // Fetch the sequence.
    //
    string  seq;
    uint8_t j;

    seq.reserve(this->aln->core.l_qseq);

    for (int i = 0; i < this->aln->core.l_qseq; i++) {
        j = bam_seqi(bam_get_seq(this->aln), i);
        switch(j) {
        case 1:
            seq += 'A';
            break;
        case 2:
            seq += 'C';
            break;
        case 4:
            seq += 'G';
            break;
        case 8:
            seq += 'T';
            break;
        case 15:
            seq += 'N';
            break;
        }
    }

    //
    // Fetch the quality score.
    //
    string   qual;
    uint8_t *q = bam_get_qual(this->aln);
    for (int i = 0; i < this->aln->core.l_qseq; i++) {
        qual += char(int(q[i]) + 33);
    }

    string chr = this->chrs[this->aln->core.tid];

    //
    // Attempt to parse the query name for this read.
    //

    s = Seq((const char *) bam_get_qname(this->aln), seq.c_str(), qual.c_str());

    return 1;
}

#endif // __BAMUNALIGNEDI_H__
