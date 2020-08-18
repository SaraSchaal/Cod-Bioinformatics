// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __CLEAN_H__
#define __CLEAN_H__

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "input.h"
#include "kmers.h"

enum fastqt   {generic_fastq, illv1_fastq, illv2_fastq};

enum barcodet {null_null,
               null_inline,   null_index,
               inline_null,   index_null,
               inline_inline, index_index,
               inline_index,  index_inline};
enum seqt {single_end, paired_end};

typedef unordered_map<string, vector<int>, std::hash<string> > AdapterHash;

extern uint     min_bc_size_1, max_bc_size_1, min_bc_size_2, max_bc_size_2;
extern int      barcode_dist_1, barcode_dist_2;
extern barcodet barcode_type;
extern uint     truncate_seq;
extern double   win_size;
extern bool     paired;
extern bool     recover;

class BarcodePair {
public:
    string se;    // Single-end barcode.
    string pe;    // Paired-end barcode.
    string name;  // Filename to open for this barcode combination.

    BarcodePair()
    {
        this->se = "";
        this->pe = "";
    }
    BarcodePair(char *p)
    {
        this->se = string(p);
        this->pe = "";
    }
    BarcodePair(char *p, char *q, char *n)
    {
        if (p != NULL)
            this->se = string(p);
        if (q != NULL)
            this->pe = string(q);
        if (n != NULL)
            this->name = string(n);
    }
    BarcodePair(string se, string pe, string name)
    {
        this->se = se;
        this->pe = pe;
        this->name = name;
    }
    BarcodePair(string se)
    {
        this->se = se;
        this->pe = "";
    }
    void set(char *p, char *q)
    {
        this->se = string(p);
        this->pe = string(q);
    }
    void set(char *p)
    {
        this->se = string(p);
        this->pe = "";
    }
    void set(string p, string q)
    {
        this->se = p;
        this->pe = q;
    }
    void set(string p)
    {
        this->se = p;
        this->pe = "";
    }
    string str()
    {
        if (this->pe.length() > 0)
            return string(this->se + "-" + this->pe);
        else
            return this->se;
    }
    bool name_exists()
    {
        if (this->name.length() > 0)
            return true;
        return false;
    }
    friend bool operator<(const BarcodePair &lhs, const BarcodePair &rhs)
    {
        if (lhs.se < rhs.se)
            return true;
        else if (lhs.se == rhs.se && lhs.pe < rhs.pe)
            return true;
        else
            return false;
    }
    friend bool operator==(const BarcodePair &lhs, const BarcodePair &rhs)
    {
        return (lhs.se == rhs.se && lhs.pe == rhs.pe);
    }
    friend ofstream& operator<<(ofstream &out, const BarcodePair &bp)
    {
        if (bp.pe.length() > 0)
            out << bp.se << "-" << bp.pe;
        else
            out << bp.se;
        return out;
    }
};

class RawRead {
public:
    fastqt fastq_type;
    char  *inline_bc;
    char  *index_bc;
    char  *se_bc;
    char  *pe_bc;
    char  *machine;
    int    run;
    int    lane;
    int    tile;
    int    x;
    int    y;
    int    index;
    int    read;
    char  *seq;
    char  *phred;
    int   *int_scores;
    bool   filter;
    int    inline_bc_len;
    int    retain;
    uint   size;
    uint   len;
    double win_len;
    double stop_pos;

    RawRead(uint buf_len, int read, int barcode_size, double win_size) {
        this->inline_bc     = new char[id_len  + 1];
        this->index_bc      = new char[id_len  + 1];
        this->machine       = new char[id_len  + 1];
        this->seq           = new char[buf_len + 1];
        this->phred         = new char[buf_len + 1];
        this->int_scores    = new  int[buf_len];
        this->size          = buf_len + 1;
        this->read          = read;

        this->retain        = 1;
        this->inline_bc_len = 0;
        this->tile          = 0;
        this->run           = 0;
        this->lane          = 0;
        this->x             = 0;
        this->y             = 0;
        this->index         = 0;
        this->len           = 0;

        this->inline_bc[0] = '\0';
        this->index_bc[0]  = '\0';
        this->machine[0]   = '\0';
        this->seq[0]       = '\0';
        this->phred[0]     = '\0';

        this->set_len(buf_len);

        this->se_bc = NULL;
        this->pe_bc = NULL;
        if (this->read == 1) {
            switch(barcode_type) {
            case index_inline:
                this->se_bc = this->index_bc;
                this->pe_bc = this->inline_bc;
                break;
            case inline_index:
                this->se_bc = this->inline_bc;
                this->pe_bc = this->index_bc;
                this->inline_bc_len = barcode_size;
                break;
            case inline_null:
            case inline_inline:
                this->se_bc = this->inline_bc;
                this->inline_bc_len = barcode_size;
                break;
            case index_null:
            case index_index:
                this->se_bc = this->index_bc;
                break;
            default:
                break;
            }
        } else if (this->read == 2) {
            switch(barcode_type) {
            case null_inline:
            case inline_inline:
            case index_inline:
                this->pe_bc = this->inline_bc;
                this->inline_bc_len = barcode_size;
                break;
            case index_index:
            case inline_index:
                this->pe_bc = this->index_bc;
                break;
            default:
                break;
            }
        }
    }
    ~RawRead() {
        delete [] this->inline_bc;
        delete [] this->index_bc;
        delete [] this->machine;
        delete [] this->seq;
        delete [] this->phred;
        delete [] this->int_scores;
    }
    int resize(int size) {
        delete [] this->seq;
        delete [] this->phred;
        delete [] this->int_scores;
        this->size  = size;
        this->seq   = new char[this->size];
        this->phred = new char[this->size];
        this->int_scores = new int[this->size - 1];

        this->set_len(size - 1);

        return 0;
    }
    int set_len(uint buf_len) {
        if (buf_len == this->len)
            return 0;

        if (buf_len > this->size - 1)
            buf_len = this->size - 1;

        this->seq[buf_len]   = '\0';
        this->phred[buf_len] = '\0';

        //
        // Set the parameters for checking read quality later in processing.
        // Window length is 15% (rounded) of the sequence length.
        //
        this->len      = buf_len - this->inline_bc_len;
        this->win_len  = round((double) this->len * win_size);

        if (this->win_len < 1)
            this->win_len = 1;

        this->len     += this->inline_bc_len;
        this->stop_pos = this->len - this->win_len;

        return 0;
    }
};

int  parse_illumina_v1(const char *);
int  parse_illumina_v2(const char *);
int  parse_input_record(Seq *, RawRead *);
int  rev_complement(char *, int, bool);
int  reverse_qual(char *, int, bool);

bool correct_barcode(set<string> &, RawRead *, seqt, int);

int  filter_adapter_seq(RawRead *, char *, int, AdapterHash &, int, int, int);
int  init_adapter_seq(int, char *, int &, AdapterHash &);

int  check_quality_scores(RawRead *, int, int, int, int);

//
// Templated function to process barcodes.
//
template<typename fhType>
int
process_barcode(RawRead *href_1, RawRead *href_2, BarcodePair &bc,
                map<BarcodePair, fhType *> &fhs,
                set<string> &se_bc, set<string> &pe_bc,
                map<BarcodePair, map<string, long> > &barcode_log, map<string, long> &counter)
{
    if (barcode_type == null_null)
        return 0;

    //
    // Is this a legitimate barcode? The barcode passed into this function is the maximally long
    // barcode. If we fail to find a match at maximum length, step down to minimum length and
    // continue to search for a match.
    //
    char *p;
    char  bc_1[id_len];
    char  bc_2[id_len];
    strcpy(bc_1, bc.se.c_str());
    strcpy(bc_2, bc.pe.c_str());

    bool valid_se_bc = false;
    bool valid_pe_bc = false;

    p = bc_1 + max_bc_size_1; // Point p at the end of string NULL.
    for (uint i = max_bc_size_1; i >= min_bc_size_1; i--)
        if (se_bc.count(bc_1) > 0) {
            valid_se_bc = true;
            break;
        } else {
            p--;
            *p = '\0';
        }
    if (pe_bc.size() > 0) {
        p = bc_2 + max_bc_size_2; // Point p at the end of string NULL.
        for (uint i = max_bc_size_2; i >= min_bc_size_2; i--)
            if (pe_bc.count(bc_2) > 0) {
                valid_pe_bc = true;
                break;
            } else {
                p--;
                *p = '\0';
            }
    }
    if (valid_se_bc == true && valid_pe_bc == true)
        bc.set(bc_1, bc_2);
    else if (valid_se_bc == true)
        bc.se = bc_1;
    else if (valid_pe_bc == true)
        bc.pe = bc_2;

    //
    // Log the barcodes we receive.
    //
    if (barcode_log.count(bc) == 0) {
        barcode_log[bc]["noradtag"] = 0;
        barcode_log[bc]["total"]    = 0;
        barcode_log[bc]["low_qual"] = 0;
        barcode_log[bc]["retained"] = 0;
    }
    barcode_log[bc]["total"] += paired ? 2 : 1;

    //
    // If we have a perfectly matching barcode, set the barcode and length in the right places.
    //
    if (pe_bc.size() > 0 && valid_se_bc == true && valid_pe_bc == true) {
        if (fhs.count(bc) > 0) {
            if (paired) {
                strcpy(href_1->se_bc, bc_1);
                strcpy(href_2->pe_bc, bc_2);
            } else {
                strcpy(href_1->se_bc, bc_1);
                strcpy(href_1->pe_bc, bc_2);
            }

            if (barcode_type == inline_index ||
                barcode_type == inline_inline)
                href_1->inline_bc_len = strlen(bc_1);
            if (barcode_type == index_inline ||
                barcode_type == inline_inline)
                href_2->inline_bc_len = strlen(bc_2);
            return 0;
        }

    } else if (valid_se_bc == true) {
        strcpy(href_1->se_bc, bc_1);
        if (barcode_type == inline_null ||
            barcode_type == inline_index ||
            barcode_type == inline_inline)
            href_1->inline_bc_len = strlen(bc_1);

    } else if (valid_pe_bc == true) {
        if (paired)
            strcpy(href_2->pe_bc, bc_2);
        else
            strcpy(href_1->pe_bc, bc_2);

        if (barcode_type == index_inline ||
            barcode_type == inline_inline)
            href_2->inline_bc_len = strlen(bc_2);
    }

    //
    // Try to correct the barcode.
    //
    BarcodePair old_barcode = bc;
    bool se_correct = false;
    bool pe_correct = false;

    if (paired) {
        if (se_bc.count(bc.se) == 0)
            se_correct = correct_barcode(se_bc, href_1, single_end, barcode_dist_1);
        if (pe_bc.size() > 0 && pe_bc.count(bc.pe) == 0)
            pe_correct = correct_barcode(pe_bc, href_2, paired_end, barcode_dist_2);

        if (se_correct)
            bc.se = string(href_1->se_bc);
        if (pe_bc.size() > 0 && pe_correct)
            bc.pe = string(href_2->pe_bc);

        //
        // After correcting the individual barcodes, check if the combination is valid.
        //
        if (fhs.count(bc) == 0) {
            counter["ambiguous"] += 2;
            href_1->retain = 0;
            href_2->retain = 0;
        }

    } else {
        if (se_bc.count(bc.se) == 0)
            se_correct = correct_barcode(se_bc, href_1, single_end, barcode_dist_1);
        if (pe_bc.size() > 0 && pe_bc.count(bc.pe) == 0)
            pe_correct = correct_barcode(pe_bc, href_1, paired_end, barcode_dist_2);

        if (se_correct)
            bc.se = string(href_1->se_bc);
        if (pe_bc.size() > 0 && pe_correct)
            bc.pe = string(href_1->pe_bc);

        if (fhs.count(bc) == 0) {
            counter["ambiguous"]++;
            href_1->retain = 0;
        }
    }

    if (href_1->retain && (se_correct || pe_correct)) {
        counter["recovered"] += paired ? 2 : 1;
        barcode_log[old_barcode]["total"] -= paired ? 2 : 1;
        if (barcode_log.count(bc) == 0) {
            barcode_log[bc]["total"]    = 0;
            barcode_log[bc]["retained"] = 0;
            barcode_log[bc]["low_qual"] = 0;
            barcode_log[bc]["noradtag"] = 0;
        }
        barcode_log[bc]["total"] += paired ? 2 : 1;
    }

    return 0;
}

#endif // __CLEAN_H__
