// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __SAMI_H__
#define __SAMI_H__

//
// Code to parse Sam format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> ... <sequence> <phred quality score> ...
//
// One record per line.
//
#include "stacks.h"
#include "input.h"
# include "aln_utils.h"
#include "BamI.h"

class Sam: public Input {

public:
    Sam(const char *path);
    ~Sam() {}
    Seq *next_seq();
    int  next_seq(Seq& s);

private:
    static const size_t n_mandatory_fields = 11;
};

Sam::Sam(const char *path)
: Input(path)
{
    while(fh.peek() == '@') {
        fh.getline(line, max_len);
    }
}

Seq* Sam::next_seq() {
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int Sam::next_seq(Seq& s) {
    vector<string> parts;
    int  flag  = 0;

    if(!fh.getline(line, max_len))
        return false;

    int len = strlen(line);
    if (line[len - 1] == '\r')
        line[len - 1] = '\0';

    parse_tsv(line, parts);
    if (parts.size() < n_mandatory_fields) {
        cerr << "Error: Malformed SAM record:\n" << line << "\n";
        throw exception();
    }

    //
    // According to SAM spec FLAGs are the second field,
    // if FLAG bit 0x4 is set, sequence is not mapped.
    //
    flag = atoi(parts[1].c_str());

    //
    // Parse the type of the record.
    //
    AlnT aln_type;
    if (flag & BAM_FUNMAP)
        aln_type = AlnT::null;
    else if (flag & BAM_FSECONDARY)
        aln_type = AlnT::secondary;
    else if (flag & BAM_FSUPPLEMENTARY)
        aln_type = AlnT::supplementary;
    else
        aln_type = AlnT::primary;

    if (aln_type == AlnT::null) {
        s = Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str()); // Read ID, Sequence, Quality
    } else {
        //
        // Check which strand this is aligned to:
        //
        strand_type strand = flag & BAM_FREVERSE ? strand_minus : strand_plus;

        //
        // Parse the alignment CIGAR string.
        //
        vector<pair<char, uint> > cigar;
        parse_cigar(parts[5].c_str(), cigar, true);
        if (strand == strand_minus)
            std::reverse(cigar.begin(), cigar.end());

        //
        // If the read was aligned on the reverse strand (and is therefore reverse complemented)
        // alter the start point of the alignment to reflect the right-side of the read, at the
        // end of the RAD cut site.
        //
        int bp = bam_find_start_bp(atoi(parts[3].c_str()), strand, cigar);
        bp--; // SAM uses 1-based genome positions.

        //
        // Calculate the percentage of the sequence that was aligned to the reference.
        //
        uint clipped = 0;
        for (auto& op : cigar)
            if (op.first == 'S')
                clipped += op.second;
        double pct_clipped = (double) clipped / parts[9].length();

        int map_qual = is_integer(parts[4].c_str());
        if (map_qual == -1)
            // Ignore malformed quality.
            map_qual = 255;

        s = Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), // Read ID, Sequence, Quality
                parts[2].c_str(), bp, strand, aln_type, pct_clipped, map_qual); // Chromosome, etc.

        if (cigar.size() > 0)
            bam_edit_gaps(cigar, s.seq);
    }
    return true;
}

#endif // __SAMI_H__
