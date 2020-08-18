#include <ctime>

#include "Vcf.h"

const VcfMeta VcfMeta::predefs::info_AD ("INFO","<ID=AD,Number=R,Type=Integer,Description=\"Total Depth for Each Allele\">");
const VcfMeta VcfMeta::predefs::info_AF ("INFO","<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
const VcfMeta VcfMeta::predefs::info_DP ("INFO","<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
const VcfMeta VcfMeta::predefs::info_NS ("INFO","<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");

const VcfMeta VcfMeta::predefs::format_AD ("FORMAT","<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">");
const VcfMeta VcfMeta::predefs::format_DP ("FORMAT","<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
const VcfMeta VcfMeta::predefs::format_HQ ("FORMAT","<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
const VcfMeta VcfMeta::predefs::format_GL ("FORMAT","<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood\">");
const VcfMeta VcfMeta::predefs::format_GQ ("FORMAT","<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
const VcfMeta VcfMeta::predefs::format_GT ("FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

const VcfMeta VcfMeta::predefs::info_loc_strand (
        "INFO",
        "<ID=loc_strand,Number=1,Type=Character,Description=\"Genomic strand the corresponding Stacks locus aligns on\">"
        );

const string VcfHeader::std_fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

void VcfHeader::add_std_meta(const string& version) {
    add_meta(VcfMeta("fileformat", version));

    time_t t;
    time(&t);
    char date[9];
    strftime(date, 9, "%Y%m%d", localtime(&t));
    add_meta(VcfMeta("fileDate", date));

    add_meta(VcfMeta("source", string("\"Stacks v") + VERSION + "\""));

    add_meta(VcfMeta::predefs::info_AD);
    add_meta(VcfMeta::predefs::info_AF);
    add_meta(VcfMeta::predefs::info_DP);
    add_meta(VcfMeta::predefs::info_NS);
    add_meta(VcfMeta::predefs::format_AD);
    add_meta(VcfMeta::predefs::format_DP);
    add_meta(VcfMeta::predefs::format_HQ);
    add_meta(VcfMeta::predefs::format_GL);
    add_meta(VcfMeta::predefs::format_GQ);
    add_meta(VcfMeta::predefs::format_GT);
}

void VcfRecord::assign(const char* rec, size_t len, const VcfHeader& header) {
    assert(rec[len-1] != '\n');

    clear();
    buffer_.assign(rec, rec+len+1);

    // Parse the record.
    size_t n_exp_fields = header.samples().empty() ? 8 : 9 + header.samples().size();
    auto wrong_n_fields = [&rec, &len, &n_exp_fields](size_t obs){
        cerr << "Error: Expected VCF record to have "
             << n_exp_fields << " fields, not " << obs << "."
             << "Offending record (" << len << " characters) is:\n"
             << rec << '\n';
        throw exception();
    };

    size_t n_fields = 1;
    char* q = NULL; // The start of the current field.
    char* p = buffer_.data(); // The start of the next field.
    auto next_field = [&p, &q, &n_fields, &wrong_n_fields]() {
        q = p;
        p = strchr(p, '\t');
        if (!p)
            wrong_n_fields(n_fields);
        ++n_fields;
        *p = '\0';
        ++p;
    };

    // chrom
    next_field();
    assert(q == buffer_.data());

    // pos
    next_field();
    pos_ = q - buffer_.data();

    // id
    next_field();

    // alleles (ref + alt)
    next_field();
    allele0_ = q - buffer_.data();
    next_field();
    while((q = strchr(q, ','))) {
        *q = '\0'; // ALT alleles become null-separated.
        ++q;
    }

    // qual
    next_field();
    qual_ = q - buffer_.data();

    // filters
    next_field();

    info0_ = p - buffer_.data();
    if (header.samples().empty()) {
        // info (last field)
        while((p = strchr(p, ';'))) {
            *p = '\0'; // INFO fields become null-separated.
            ++p;
        }
    } else {
        // info
        next_field();
        while((q = strchr(q, ';'))) {
            *q = '\0';
            ++q;
        }

        // format
        next_field();
        format0_ = q - buffer_.data();
        while((q = strchr(q, ':'))) {
            *q = '\0'; // FORMAT fields become null-separated.
            ++q;
        }

        // samples (last fields)
        sample0_ = p - buffer_.data();
        while((p = strchr(p, '\t'))) {
            *p = '\0';
            ++p;
            ++n_fields;
        }
    }
    if (n_fields != n_exp_fields)
        wrong_n_fields(n_fields);
}

string VcfRecord::util::fmt_info_af(const vector<double>& alt_freqs) {
    stringstream ss;
    ss << std::setprecision(3);
    join(alt_freqs, ',', ss);
    return string("AF=") + ss.str();
}

string VcfRecord::util::fmt_gt_gl(const vector<Nt2>& alleles, const GtLiks& liks) {
    stringstream ss;
    ss << std::fixed << std::setprecision(2);

    assert(!alleles.empty());
    vector<double> v;
    v.reserve(n_possible_genotypes(alleles.size()));
    for (size_t a2=0; a2<alleles.size(); ++a2) {
        Nt2 a2nt = alleles[a2];
        for (size_t a1=0; a1<=a2; ++a1) {
            Nt2 a1nt = alleles[a1];
            v.push_back(liks.at(a1nt, a2nt) / log(10));
        }
    }
    join(v, ',', ss);
    return ss.str();
}

size_t VcfRecord::util::parse_gt_dp(
    const char* gt_str,
    size_t dp_index
) {
    const char* p = VcfRecord::util::find_gt_subfield(gt_str, dp_index);
    if (p == NULL || *p == '.') {
        cerr << "Error: DP field is missing.\n";
        throw exception();
    }
    char* tmp;
    long dp = strtol(p, &tmp, 10);
    if (tmp == p || dp < 0) {
        cerr << "Error: Bad DP field.\n";
        throw exception();
    }
    return dp;
}

Counts<Nt2>
VcfRecord::util::parse_gt_ad(
    const char* gt_str,
    size_t ad_index,
    const vector<Nt2>& alleles
) {
    Counts<Nt2> depths;
    const char* p = VcfRecord::util::find_gt_subfield(gt_str, ad_index);
    if (p == NULL || *p == '.') {
        cerr << "Error: AD field is missing.\n";
        throw exception();
    }
try {
    size_t i = 0;
    char* end;
    --p;
    do {
        ++p;
        long ad = strtol(p, &end, 10);
        if (end == p || ad < 0)
            throw exception();
        depths.increment(alleles.at(i), ad);
        ++i;
        p = end;
    } while (*p == ',');
    if (i != alleles.size()) {
        cerr << "Error: Wrong number of values (expected "
             << alleles.size() << ", found" << i << ").\n";
        throw exception();
    }
    return depths;
} catch (exception&) {
    cerr << "Error: Bad AD field in '" << gt_str << "'.\n";
    throw;
}}

uint8_t VcfRecord::util::parse_gt_gq(
    const char* gt_str,
    size_t gq_index
) {
    const char* p = VcfRecord::util::find_gt_subfield(gt_str, gq_index);
    if (p == NULL || *p == '.') {
        cerr << "Error: GQ field is missing.\n";
        throw exception();
    }
    char* tmp;
    long gq = strtol(p, &tmp, 10);
    if (tmp == p || gq < 0 || gq > UINT8_MAX) {
        cerr << "Error: Bad GQ field.\n";
        throw exception();
    }
    return gq;
}

GtLiks VcfRecord::util::parse_gt_gl(
    const char* gt_str,
    size_t gl_index,
    const vector<Nt2>& alleles
) {
    GtLiks liks;
    const char* p = VcfRecord::util::find_gt_subfield(gt_str, gl_index);
    if (p == NULL) {
        cerr << "Error: GL field is missing.\n";
        throw exception();
    }
try {
    size_t n_values_found = 0;
    char* end;
    for (size_t a2=0; a2<alleles.size(); ++a2) {
        Nt2 a2nt = alleles[a2];
        for (size_t a1=0; a1<=a2; ++a1) {
            Nt2 a1nt = alleles[a1];
            double lik = strtod(p, &end);
            if (end == p)
                throw exception();
            lik *= log(10); // Back to base e.
            liks.set(a1nt, a2nt, lik);
            // Move to the next field.
            ++n_values_found;
            p = end;
            if (*p == ',') {
                ++p;
            } else if (*p != ':' && *p != '\0') {
                throw exception();
            } else if (n_values_found != n_possible_genotypes(alleles.size())) {
                cerr << "Error: Wrong number of values (expected "
                     << n_possible_genotypes(alleles.size()) << ", found "
                     << n_values_found << ").\n";
                throw exception();
            }
        }
    }
    return liks;
} catch (exception&) {
    cerr << "Error: VCF: Illegal GL format: "
        << (gl_index+1) << "th field in '" << gt_str << "'.\n";
    throw;
}}

VcfParser::VcfParser(const string& path) : file_(path), header_() {
    FileT ftype = guess_file_type(path);
    if (ftype != FileT::vcf && ftype != FileT::gzvcf) {
        cerr << "Error: File '" << path << "' : expected '.vcf(.gz)' suffix.\n";
        throw exception();
    }
    read_header();
}

int
VcfParser::open(string &path) {
    this->file_.open(path);

    FileT ftype = guess_file_type(path);

    if (ftype != FileT::vcf && ftype != FileT::gzvcf) {
        cerr << "Error: File '" << path << "' : expected '.vcf(.gz)' suffix.\n";
        throw exception();
    }

    read_header();

    return 0;
}

void
VcfParser::read_header()
{
    const char* line = NULL;
    size_t len;

    auto malformed = [this, &line] () {
        cerr << "Error: Malformed VCF header."
             << " At line " << file_.line_number() << " in file '" << file_.path() << "'";
        if (line)
            cerr << ":\n" << line << "\n";
        else
            cerr << ".\n";
        throw exception();
    };

    // Meta header lines.
    bool not_eof;
    while((not_eof = file_.getline(line, len)) && strncmp(line, "##", 2) == 0) {
        const char* equal = strchr(line+2, '=');
        if(!equal) {
            // Skip header lines missing the '='.
            continue;
        }
        header_.add_meta(VcfMeta(string(line+2, equal), string(equal+1)));
    }
    if (!not_eof)
        malformed();

    // Final header line.
    if(strncmp(line, VcfHeader::std_fields.c_str(), VcfHeader::std_fields.length()) != 0)
        malformed();

    // Parse sample names, if any.
    if(len > VcfHeader::std_fields.length()) {
        const char* p = line + VcfHeader::std_fields.length();
        const char* end = line + len;

        // Check that FORMAT is present.
        const char format[] = "\tFORMAT";
        const size_t format_len = sizeof(format) - 1;
        if(strncmp(p, format, format_len) != 0)
            malformed();
        p += format_len;
        if (*p != '\t')
            malformed();

        do {
            ++p;
            const char* next = strchr(p, '\t');
            if (!next)
                next = end;
            header_.add_sample(string(p, next));
            p = next;
        } while (p != end);
    }
}
