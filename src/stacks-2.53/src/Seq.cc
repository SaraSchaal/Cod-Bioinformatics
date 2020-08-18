#include <utility>

#include "constants.h"
#include "utils.h"
#include "Seq.h"

const char PhyLoc::empty_str[1] = {'\0'};

PhyLoc::PhyLoc(const string& s) : PhyLoc() {
    try {
        // Chromosome.
        const char* p = strchr(s.c_str(), ':');
        if (p == NULL || p == s.c_str())
            throw exception();
        chr_ = new char[p-s.c_str()+1];
        strncpy(chr_, s.c_str(), p-s.c_str());
        chr_[p-s.c_str()] = '\0';

        // BP.
        ++p;
        char* end;
        bp = strtol(p, &end, 10) - 1;
        if (end == p)
            throw exception();

        // Strand.
        if (*end != '\0') {
            if (*end != ':')
                throw exception();
            p = end + 1;
            if (!( (*p == '+' || *p == '-') && *(p+1) == '\0' ) )
                throw exception();
            strand = (*p == '+') ? strand_plus : strand_minus;
        }
    } catch (exception&) {
        cerr << "Error: Malformed genomic position '" << s << "'.\n";
        throw;
    }
}

Seq::Seq()
    : id (NULL)
    , capacity (-1)
    , seq (NULL)
    , qual (NULL)
    , aln_type (AlnT::null)
    , pct_clipped (0.0)
    , map_qual (255)
    , loc_str  (NULL)
{}

Seq::Seq(const Seq& other)
    : capacity(other.capacity)
    , comment(other.comment)
    , loc(other.loc)
{
    if (other.id != NULL) {
        id = new char[strlen(other.id)+1];
        strcpy(id, other.id);
    } else {
        id = NULL;
    }
    if (capacity < 0) {
        capacity = std::max(
            (other.seq  != NULL ? strlen(other.seq)  : -1),
            (other.qual != NULL ? strlen(other.qual) : -1)
            );
    }
    if (other.seq != NULL) {
        seq = new char[capacity + 1];
        assert(long(strlen(other.seq)) <= capacity);
        strcpy(seq, other.seq);
    } else {
        seq = NULL;
    }
    if (other.qual != NULL) {
        qual = new char[capacity + 1];
        assert(long(strlen(other.qual)) <= capacity);
        strcpy(qual, other.qual);
    } else {
        qual = NULL;
    }
    if (other.loc_str != NULL) {
        loc_str = new char[strlen(other.loc_str)+1];
        strcpy(loc_str, other.loc_str);
    } else {
        loc_str = NULL;
    }

    pct_clipped  = other.pct_clipped;
    aln_type = other.aln_type;
    map_qual = other.map_qual;
}

Seq::Seq(const char *id, const char *seq)
    : Seq()
{
    this->id = new char[strlen(id) + 1];
    strcpy(this->id, id);

    capacity = strlen(seq);
    this->seq = new char[capacity + 1];
    strcpy(this->seq, seq);
}

Seq::Seq(const char *id, const char *seq, const char *qual)
    : Seq(id, seq)
{
    this->qual = new char[capacity + 1];
    // if (strlen(qual) > capacity)
    //     throw std::invalid_argument("Seq::Seq(id, seq, qual)");
    strncpy(this->qual, qual, capacity); // (trucates)
    this->qual[capacity] = '\0';
}

Seq::Seq(const char *id, const char *seq, const char *qual,
         const char *chr, uint bp, strand_type strand
)
    : Seq(id, seq, qual)
{
    this->loc.set(chr, bp, strand);
    this->aln_type = AlnT::primary;

    // Reverse complement sequences from the negative strand
    if (strand == strand_minus) {
        rev_comp_inplace(this->seq);
        if (this->qual != NULL)
            std::reverse(this->qual, this->qual + strlen(this->qual));
    }

    // Set `loc_str`.
    this->loc_str = new char[strlen(chr)  + 15];
    sprintf(this->loc_str, "%s|%d|%c", chr, bp, strand == strand_plus ? '+' : '-');
}

Seq::Seq(const char *id, const char *seq, const char *qual,
         const char *chr, uint bp, strand_type strand,
         AlnT aln_type, double pct_clipped, int map_qual
)
    : Seq(id, seq, qual, chr, bp, strand)
{
    this->aln_type = aln_type;
    this->pct_clipped  = pct_clipped;
    this->map_qual = map_qual;
}

void swap(Seq& s1, Seq& s2) {
    std::swap(s1.id, s2.id);
    std::swap(s1.capacity, s2.capacity);
    std::swap(s1.seq, s2.seq);
    std::swap(s1.qual, s2.qual);
    s1.comment.swap(s2.comment);
    std::swap(s1.loc_str, s2.loc_str);
    std::swap(s1.aln_type, s2.aln_type);
    std::swap(s1.pct_clipped, s2.pct_clipped);
    swap(s1.loc, s2.loc); // (ADL)
    std::swap(s1.map_qual, s2.map_qual);
}
