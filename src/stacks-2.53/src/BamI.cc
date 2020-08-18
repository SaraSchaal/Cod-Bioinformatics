#include <cstring>
#include <string>
#include <map>

#include "constants.h"
#include "BamI.h"

using namespace std;

const uint32_t o = UINT32_MAX;
const uint32_t cigar_c2i[128] = {
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,7,o,o,

    o,o,9,o, 2,o,o,o, 5,1,o,o, o,0,3,o,
    6,o,o,4, o,o,o,o, 8,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
};

void BamRecord::assign(
        const string& name,
        uint16_t flg,
        int32_t chr_index,
        int32_t aln_pos,
        const vector<pair<char,uint>>& cig,
        const DNASeq4& seq,
        size_t read_group
        ) {
    if (empty())
        reinit();

    // bam1_t::core
    r_->core.tid = chr_index;
    r_->core.pos = aln_pos;
    r_->core.bin = 0; // No idea
    r_->core.qual = 255;
    r_->core.l_qname = name.length() + 1; // `l_qname` includes the trailing \0.
    r_->core.flag = flg;
    r_->core.n_cigar = cig.size();
    r_->core.l_qseq = seq.length();
    r_->core.mtid = -1;
    r_->core.mpos = -1;
    r_->core.isize = -1; // No idea

    // bam1_t::data
    // Htslib says: "bam1_t::data -- all variable-length data, concatenated;
    // structure: qname-cigar-seq-qual-aux, concatenated".

    // Prepare the `aux` data.
    string rg = string() + "RG" + "Z" + to_string(read_group);

    // Determine the length of `data`.
    size_t l_aux = rg.length() + 1;
    r_->l_data = r_->core.l_qname + r_->core.n_cigar*sizeof(uint32_t) + seq.nbytes() + seq.length() + l_aux;
    if ((uint)r_->l_data > r_->m_data) {
        if (r_->data != NULL)
            free(r_->data);
        r_->m_data = r_->l_data;
        r_->data = (uchar*) malloc(r_->m_data);
    }

    // Fill the data array.
    uchar* p = r_->data;
    //qname
    strcpy((char*)p, name.c_str());
    p += r_->core.l_qname;
    //cigar
    for (const pair<char, uint>& op : cig) {
        // Cigars are uint32_t's with the length on the 28 high bits & op on the low 4 bits.
        *(uint32_t*)p = (uint32_t(op.second) <<BAM_CIGAR_SHIFT) | cigar_c2i[size_t(op.first)];
        p += sizeof(uint32_t);
    }
    //seq & qual
    memcpy(p, seq.vdata(), seq.nbytes());
    p += seq.nbytes();
    memset(p, 0xFF, seq.length());
    p += seq.length();
    //aux
    memcpy(p, rg.c_str(), rg.length()+1);

    // bam1_t::core.bin
    // c.f. `sam_parse1()`; I have no idea what this is.
    uint32_t* cigar = (uint32_t*)(r_->data + r_->core.l_qname);
    r_->core.bin = hts_reg2bin(r_->core.pos, r_->core.pos + bam_cigar2rlen(r_->core.n_cigar, cigar), 14, 5);
}

BamHeader::BamHeader(const char* text, size_t len): h_() {
    h_ = sam_hdr_parse(len+1, text); // null-terminated
    if (h_ == NULL)
        throw ios::failure("sam_hdr_parse");
    h_->l_text = len+1;
    h_->text = (char*) malloc(h_->l_text);
    strcpy(h_->text, text);
    init_sdict();
}

void BamHeader::init_sdict() {
    if(h_->n_targets == 0)
        DOES_NOT_HAPPEN;
    int32_t i = bam_name2id(h_, h_->target_name[0]);
    assert(h_->sdict);
    assert(i == 0);
}

void BamHeader::check_same_ref_chroms(
        const BamHeader& h1,
        const BamHeader& h2
) {
    if (h1.n_ref_chroms() != h2.n_ref_chroms()) {
        cerr << "Error: Headers have different number of chromosomes.\n";
        throw exception();
    }
    for (size_t i=0; i<h1.n_ref_chroms(); ++i) {
        if (strcmp(h1.chrom_str(i), h2.chrom_str(i)) != 0) {
            cerr << "Error: Conflicting names for the " << i+1 << "th chromosome, '"
                 << h1.chrom_str(i) << "' and '" << h2.chrom_str(i) << "'.\n";
            throw exception();
        }
        if (h1.chrom_len(i) != h2.chrom_len(i)) {
            cerr << "Error: " << i+1 << "th chromosome has lengths "
                 << h1.chrom_len(i) << " and " << h2.chrom_len(i) << ".\n";
            throw exception();
        }
    }
}

Bam::Bam(const char *path)
:
    Input(),
    bam_fh(hts_open(path, "r")),
    hdr(),
    eof_(false),
    n_records_read_(0),
    prev_chrom_(0),
    prev_pos_(0)
{
    this->path   = string(path);
    check_open(bam_fh, path);
    if (bam_fh->format.format != bam) {
        cerr << "Error: '" << path << "':";
        if (bam_fh->format.format == sam)
            cerr << " this is a SAM file (and BAM was specified).\n";
        else
            cerr << " not a BAM file.\n";
        throw exception();
    }
    hdr.reinit(bam_fh);
};

Bam::Bam(const string& path, BamHeader&& header)
:
    Input(),
    bam_fh(hts_open(path.c_str(), "wb")),
    hdr(move(header)),
    eof_(false),
    n_records_read_(0),
    prev_chrom_(0),
    prev_pos_(0)
{
    this->path   = path;
    check_open(bam_fh, path);

    // Write the header.
    int rv = bam_hdr_write(bam_fh->fp.bgzf, hdr.hts());
    if (rv != 0) {
        cerr << "Error: Writing of BAM header failed (`bam_hdr_write()`returned " << rv << ").\n";
        throw ios::failure("bam_hdr_write");
    }
}

void Bam::check_open(const htsFile* bam_f, const string& path) {
    if (bam_f == NULL) {
        #pragma omp critical (bam_check_open)
        {
            cerr << "Error: Failed to open BAM file '" << path << "'.\n";
        }
        throw exception();
    }
}

Seq *
Bam::next_seq()
{
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int
Bam::next_seq(Seq& s)
{
    //
    // Read a record
    //
    BamRecord rec;
    if (!next_record(rec))
        return false;

    //
    // Fetch the sequence.
    //
    string  seq;
    seq.reserve(rec.hts()->core.l_qseq);
    for (int i = 0; i < rec.hts()->core.l_qseq; i++) {
        uint8_t j = bam_seqi(bam_get_seq(rec.hts()), i);
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
        default:
            DOES_NOT_HAPPEN;
            break;
        }
    }

    //
    // Fetch the quality score.
    //
    string   qual;
    uint8_t *q = bam_get_qual(rec.hts());
    for (int i = 0; i < rec.hts()->core.l_qseq; i++) {
        qual += char(int(q[i]) + 33);
    }

    AlnT aln_type;
    if (rec.is_unmapped())
        aln_type = AlnT::null;
    else if (rec.is_secondary())
        aln_type = AlnT::secondary;
    else if (rec.is_supplementary())
        aln_type = AlnT::supplementary;
    else
        aln_type = AlnT::primary;

    if (aln_type == AlnT::null) {
        s = Seq(rec.qname(), seq.c_str(), qual.c_str());
    } else {
        //
        // Check which strand this is aligned to:
        //   SAM reference: FLAG bit 0x10 - sequence is reverse complemented
        //
        strand_type strand = rec.is_rev_compl() ? strand_minus : strand_plus;

        //
        // Parse the alignment CIGAR string.
        // If aligned to the negative strand, sequence has been reverse complemented and
        // CIGAR string should be interpreted in reverse.
        //
        Cigar cigar = rec.cigar();
        if (strand == strand_minus)
            std::reverse(cigar.begin(), cigar.end());

        //
        // If the read was aligned on the reverse strand (and is therefore reverse complemented)
        // alter the start point of the alignment to reflect the right-side of the read, at the
        // end of the RAD cut site.
        //
        uint bp = bam_find_start_bp(rec.pos(), strand, cigar);

        //
        // Calculate the percentage of the sequence that was aligned to the reference.
        //
        uint clipped = 0;
        for (auto& op : cigar)
            if (op.first == 'S')
                clipped += op.second;
        double pct_clipped = (double) clipped / seq.length();

        string name = rec.qname();
        if (rec.is_read1())
            name += "/1";
        else if (rec.is_read2())
            name += "/2";
        s = Seq(name.c_str(), seq.c_str(), qual.c_str(),
                hdr.chrom_str(rec.chrom()), bp, strand,
                aln_type, pct_clipped, rec.mapq());

        if (cigar.size() > 0)
            bam_edit_gaps(cigar, s.seq);
    }

    return true;
}

int
bam_find_start_bp(int aln_bp, strand_type strand, const Cigar& cigar)
{
    if (strand == strand_plus) {
        if (cigar.at(0).first == 'S')
            aln_bp -= cigar.at(0).second;
    } else {
        // assert(strand == strand_minus);
        for (uint i = 0; i < cigar.size(); i++)  {
            char op   = cigar[i].first;
            uint dist = cigar[i].second;

            switch(op) {
            case 'I':
            case 'H':
                break;
            case 'S':
                if (i < cigar.size() - 1)
                    aln_bp += dist;
                break;
            case 'M':
            case '=':
            case 'X':
            case 'D':
            case 'N':
                aln_bp += dist;
                break;
            default:
                break;
            }
        }
        aln_bp -= 1;
    }

    return aln_bp;
}

int
bam_edit_gaps(const Cigar& cigar, char *seq)
{
    char *buf;
    uint  size = cigar.size();
    char  op;
    uint  dist, bp, len, buf_len, buf_size, j, k, stop;

    len = strlen(seq);
    bp  = 0;

    buf      = new char[len + 1];
    buf_size = len + 1;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }
            break;
        case 'D':
            //
            // A deletion has occured in the read relative to the reference genome.
            // Pad the read with sufficent Ns to match the deletion, shifting the existing
            // sequence down. Trim the final length to keep the read length consistent.
            //
            k = bp >= len ? len : bp;

            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len         = strlen(buf);

            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }
            break;
        case 'I':
            //
            // An insertion has occurred in the read relative to the reference genome. Delete the
            // inserted bases and pad the end of the read with Ns.
            //
            if (bp >= len) break;

            k = bp + dist > len ? len : bp + dist;
            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len           = strlen(buf);

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }

            stop = j + dist;
            stop = stop > len ? len : stop;
            while (j < stop) {
                seq[j] = 'N';
                j++;
            }
            break;
        case 'M':
        case '=':
        case 'X':
            bp += dist;
            break;
        default:
            break;
        }
    }

    delete [] buf;

    return 0;
}
