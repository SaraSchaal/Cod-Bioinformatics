// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2019, Julian Catchen <jcatchen@illinois.edu>
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
#ifndef EXPORT_FORMATS_H
#define EXPORT_FORMATS_H

#include <iostream>
#include <utility>
#include <map>
#include <typeinfo>
#include <typeindex>

#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "ordered.h" // for "snp"
#include "mapping_utils.h"
#include "genotype_dictionaries.h"

void tally_complete_haplotypes(
        Datum const*const* data,
        size_t n_samples,
        strand_type loc_strand,
        vector<pair<const char*, size_t>>& haps_sorted_decr_freq,
        map<const char*, size_t, LessCStrs>& hap_indexes_map
        );

class Export {
 protected:
    string     _path;
    ofstream   _fh;

 public:
    Export() {}
    virtual ~Export() {}
    virtual int  open(const MetaPopInfo *) = 0;
    virtual int  write_header()    = 0;
    virtual int  write_batch(const vector<LocBin *> &) = 0;
    virtual int  post_processing() {return 0;}
    virtual void close()           {this->_fh.close();}

    bool is_hap_export();
    string tmp_path() const {return this->_path + ".part";}
    static int transpose(ifstream &ifh, vector<string> &transposed);
};

class OrderableExport : public Export {
 public:
    OrderableExport() {}
    virtual ~OrderableExport() {}
    int write_batch(const vector<LocBin*>& loci);

 protected:
    virtual int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index) = 0;
};

class GenPos {
 public:
    uint     id;
    uint     bp;
    uint     snp_index;
    loc_type type;

    GenPos(int id, int snp_index, int bp) {
        this->id        = id;
        this->snp_index = snp_index;
        this->bp        = bp;
        this->type      = snp;
    }
    GenPos(int id, int snp_index, int bp, loc_type type) {
        this->id        = id;
        this->snp_index = snp_index;
        this->bp        = bp;
        this->type      = type;
    }

    bool operator<(const GenPos& other) const {return bp < other.bp;}
};

class MarkersExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;

 public:
    MarkersExport() : _mpopi(NULL) {}
    ~MarkersExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class RawHaplotypesExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;

 public:
    RawHaplotypesExport() : _mpopi(NULL) {}
    ~RawHaplotypesExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class SumstatsExport: public Export {
    //
    // Output the locus-level summary statistics.
    //
    const MetaPopInfo *_mpopi;
    uint  _pop_cnt;

 public:
    SumstatsExport() : _mpopi(NULL), _pop_cnt(UINT_MAX) {}
    ~SumstatsExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class HapstatsExport: public Export {
    //
    // Output the locus-level haplotype statistics.
    //
    const MetaPopInfo *_mpopi;
    uint  _pop_cnt;

 public:
    HapstatsExport() : _mpopi(NULL), _pop_cnt(UINT_MAX) {}
    ~HapstatsExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class SnpDivergenceExport: public Export {
    //
    // Output the SNP-level divergence statistics.
    //
    const MetaPopInfo *_mpopi;
    vector<ofstream *> _fhs;
    OPopPair<PopPair> *_order;

 public:
    SnpDivergenceExport(ofstream &log_fh);
    ~SnpDivergenceExport() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            delete this->_fhs[i];
        delete this->_order;
    }
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &) { return 0; }
    int  write_batch_pairwise(const vector<LocBin *> &, const vector<vector<PopPair **>> &);
    void close() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            this->_fhs[i]->close();
        return;
    }

private:
    int write_site(ofstream *fh, const PopPair *pp, string chr);
};

class HapDivergenceExport: public Export {
    //
    // Output the SNP-level divergence statistics.
    //
    const MetaPopInfo *_mpopi;
    vector<ofstream *> _fhs;
    ofstream *_metapop_fh;

 public:
    HapDivergenceExport() : _mpopi(NULL), _metapop_fh(NULL) {}
    ~HapDivergenceExport() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            delete this->_fhs[i];
    }
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &) { return 0; }
    int  write_batch_pairwise(const vector<LocBin *> &, const vector<vector<HapStat *>> &, const vector<HapStat *> &);
    void close() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            this->_fhs[i]->close();
        return;
    }
};

class PlinkExport: public OrderableExport {
    const MetaPopInfo *_mpopi;
    ofstream _tmpfh;
    ifstream _intmpfh;
    string   _markers_path;
    ofstream _markers_fh;

 public:
    PlinkExport() : _mpopi(NULL) {}
    ~PlinkExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &loci);
    int  post_processing();
    void close();

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class GenePopExport: public OrderableExport {
    const MetaPopInfo *_mpopi;
    ofstream _tmpfh;

 public:
    GenePopExport() : _mpopi(NULL) {}
    ~GenePopExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  post_processing();
    void close() {this->_fh.close(); remove(this->tmp_path().c_str());}

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class GenePopHapsExport: public Export {
    const MetaPopInfo *_mpopi;
    ofstream _tmpfh;
    size_t _n_digits;

 public:
    GenePopHapsExport() : _mpopi(NULL), _n_digits(2) {}
    ~GenePopHapsExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin*>& loci);
    int  post_processing();
    void close() {this->_fh.close(); remove(this->tmp_path().c_str());}

    void set_digits(size_t n) { assert(n==2 || n==3); _n_digits=n; }

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class StructureExport: public OrderableExport {
    const MetaPopInfo *_mpopi;
    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;

 public:
    StructureExport() : _mpopi(NULL) {}
    ~StructureExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  post_processing();
    void close();

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class FineRADStructureExport: public Export {
    const MetaPopInfo *_mpopi;
    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;

 public:
    FineRADStructureExport() : _mpopi(NULL) {}
    ~FineRADStructureExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class PhylipExport: public OrderableExport {
protected:
    const MetaPopInfo *_mpopi;
    string   _log_path;
    ofstream _logfh;
    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;
    size_t   _site_index;

 public:
    PhylipExport() : _mpopi(NULL), _site_index(0) {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  post_processing();
    void close();
};

class PhylipVarExport: public PhylipExport {
 public:
    PhylipVarExport() {}
 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class PhylipFixedExport: public PhylipExport {
 public:
    PhylipFixedExport() {}
 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class PhylipVarAllExport: public Export {
    const MetaPopInfo *_mpopi;
    string   _log_path;
    ofstream _logfh;
    ofstream _parfh;
    string   _partition_path;

    size_t _seq_len;
    size_t _site_index;
    size_t _loc_cnt;
    map<int, string> _outstrs;
    
public:
    PhylipVarAllExport() : _mpopi(NULL), _seq_len(0), _site_index(1), _loc_cnt(1) {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing();
    void close();
};

class HzarExport: public OrderableExport {
    const MetaPopInfo *_mpopi;
    ofstream _tmpfh;
    string   _tmp_path;
    ifstream _intmpfh;

 public:
    HzarExport() : _mpopi(NULL) {}
    ~HzarExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header() { return 0; };
    int  post_processing();
    void close();

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class FastaLociExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;

 public:
    FastaLociExport() : _mpopi(NULL) {}
    ~FastaLociExport() {}
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
};

class FastaRawExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;

 public:
    FastaRawExport() : _mpopi(NULL) {}
    ~FastaRawExport() {}
    int open(const MetaPopInfo *mpopi);
    int write_header();
    int write_batch(const vector<LocBin *> &);
};

class FastaSamplesExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;

 public:
    FastaSamplesExport() : _mpopi(NULL) {}
    ~FastaSamplesExport() {}
    int open(const MetaPopInfo *mpopi);
    int write_header();
    int write_batch(const vector<LocBin *> &);
};

class VcfExport: public OrderableExport {
    const MetaPopInfo*_mpopi;
    VcfWriter* _writer;

 public:
    VcfExport() : _mpopi(NULL), _writer(NULL) {}
    ~VcfExport() { delete this->_writer; }
    int open(const MetaPopInfo *mpopi);
    int write_header() { return 0; }

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class VcfHapsExport: public Export {
    const MetaPopInfo*_mpopi;
    VcfWriter* _writer;

 public:
    VcfHapsExport() : _mpopi(NULL), _writer(NULL) {}
    ~VcfHapsExport() { delete this->_writer; }
    int open(const MetaPopInfo *mpopi);
    int write_batch(const vector<LocBin*>& loci);
    int write_header() { return 0; }
};

class TreemixExport: public OrderableExport {
    const MetaPopInfo*_mpopi;
    ofstream _writer;

 public:
    TreemixExport() : _mpopi(NULL), _writer(NULL) {}
    int open(const MetaPopInfo *mpopi);
    int write_header() { return 0; }

 private:
    int write_site(const CSLocus* cloc, const LocPopSum* psum, Datum const*const* datums, size_t col, size_t index);
};

class JoinMapExport: public Export {
    const MetaPopInfo*_mpopi;
    const MappingGenotypesProcessor *_gproc;
    size_t _progeny_index;

    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;
    
    map<string, string>              _marker_map;
    map<string, map<string, string>> _genotype_map;

public:
    JoinMapExport(MappingGenotypesProcessor *g) : _mpopi(NULL), _gproc(g) {}
    JoinMapExport() : _mpopi(NULL), _gproc(NULL) {}
    int    open(const MetaPopInfo *mpopi);
    int    write_batch(const vector<LocBin*>& loci);
    int    write_header() { return 0; }
    int    post_processing();
    void   close();
};

class OneMapExport: public Export {
    const MetaPopInfo*_mpopi;
    const MappingGenotypesProcessor *_gproc;
    size_t _progeny_index;

    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;

    map<string, string> _onemap_marker_types =
        {{"abxoo", "D1.11"},
         {"ooxab", "D2.16"},
         {"abxaa", "D1.10"},
         {"aaxab", "D2.15"},
         {"abxab", "B3.7"},
         {"abxac", "A.2"},
         {"abxcd", "A.1"}};

    map<string, string>              _marker_map;
    map<string, map<string, string>> _genotype_map;

 public:
    OneMapExport(MappingGenotypesProcessor *g) : _mpopi(NULL), _gproc(g) {}
    OneMapExport() : _mpopi(NULL) {}
    int    open(const MetaPopInfo *mpopi);
    int    write_batch(const vector<LocBin*>& loci);
    int    write_header() { return 0; }
    int    post_processing();
    void   close();
};

class rQTLExport: public Export {
    const MetaPopInfo*_mpopi;
    const MappingGenotypesProcessor *_gproc;
    size_t _progeny_index;

    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;

    map<string, string>              _marker_map;
    map<string, map<string, string>> _genotype_map;

 public:
    rQTLExport(MappingGenotypesProcessor *g) : _mpopi(NULL), _gproc(g) {}
    rQTLExport() : _mpopi(NULL) {}
    int    open(const MetaPopInfo *mpopi);
    int    write_batch(const vector<LocBin*>& loci);
    int    write_header();
    int    post_processing();
    void   close();
};

class MapMakerExport: public Export {
    const MetaPopInfo*_mpopi;
    const MappingGenotypesProcessor *_gproc;
    size_t _progeny_index;

    string   _tmp_path;
    ofstream _tmpfh;
    ifstream _intmpfh;

    map<string, string>              _marker_map;
    map<string, map<string, string>> _genotype_map;

 public:
    MapMakerExport(MappingGenotypesProcessor *g) : _mpopi(NULL), _gproc(g) {}
    MapMakerExport() : _mpopi(NULL) {}
    int    open(const MetaPopInfo *mpopi);
    int    write_batch(const vector<LocBin*>& loci);
    int    write_header() { return 0; }
    int    post_processing();
    void   close();
};

/*
int write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, bool);
int write_phase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_fastphase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
*/

int find_datum_allele_depths(const Datum*, int, char, char, int &, int &);
int tally_observed_haplotypes(const vector<char *> &, int, char &, char &);

#endif // EXPORT_FORMATS_H
