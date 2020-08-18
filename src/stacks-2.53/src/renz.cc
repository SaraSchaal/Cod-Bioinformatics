#include "renz.h"

using namespace std;

//
// First line of static array contains each enzyme's cut sites. Second
// line is the reverse complement of each cut site.
//
const char *aciI[]    = {"CGC", "CGG",        // C/CGC, AciI
                         "GCG", "CCG"};
const char *ageI[]    = {"CCGGT",             // A/CCGGT, AgeI
                         "ACCGG"};
const char *aluI[]    = {"CT",                // AG/CT, AluI
                         "AG"};
const char *apaLI[]   = {"TGCAC",             // G/TGCAC, ApaLI
                         "GTGCA"};
const char *apeKI[]   = {"CAGC", "CTGC",      // G/CWGC, ApeKI; W=A or T
                         "GTCG", "GACG"};
const char *apoI[]    = {"AATTC", "AATTT",    // R/AATTY, ApoI  (also known as XapI)
                         "GAATT", "AAATT"};
const char *aseI[]    = {"TAAT",              // AT/TAAT, AseI
                         "ATTA"};
const char *bamHI[]   = {"GATCC",             // G/GATCC, BamHI
                         "GGATC"};
const char *bbvCI[]   = {"TCAGC", "GGAGT",    // CC/TCAGC and GGAGT/CG, BbvCI
                         "GCTGA", "ACTCC"};
const char *bfaI[]    = {"TAG",               // C/TAG, BfaI
                         "CTA"};
const char *bfuCI[]   = {"GATC",              // /GATC, BfuCI
                         "GATC"};
const char *bgIII[]   = {"GATCT",             // A/GATCT, BgIII
                         "AGATC"};
const char *bsaHI[]   = {"CGCC", "CGTC",      // GR/CGYC, BsaHI
                         "GGCG", "GACG"};
const char *bspDI[]   = {"CGAT",              // AT/CGAT, BspDI
                         "ATCG"};
const char *bstYI[]   = {"GATCC", "GATCT",    // R/GATCY, BstYI (also known as PsuI)
                         "GGATC", "AGATC"};
const char *cac8I[]    = {"AGC", "CGC", "GGC", "TGC", // GCN/NGC, Cac8I
                          "GCT", "GCG", "GCC", "GCA"};
const char *claI[]    = {"CGAT",              // AT/CGAT, ClaI
                         "ATCG"};
const char *csp6I[]   = {"TAC",               // G/TAC, Csp6I
                         "GTA"};
const char *ddeI[]    = {"TAAG", "TCAG", "TGAG", "TTAG", // C/TNAG, DdeI
                         "CTTA", "CTGA", "CTCA", "CTAA"};
const char *dpnII[]   = {"GATC",              // GATC, DpnII
                         "GATC"};
const char *eaeI[]    = {"GGCCA", "GGCCG",    // Y/GGCCR, EaeI
                         "TGGCC", "CGGCC"};
const char *ecoRI[]   = {"AATTC",             // G/AATTC, EcoRI
                         "GAATT"};
const char *ecoRV[]   = {"ATC",               // GAT/ATC, EcoRV
                         "GAT"};
const char *ecoT22I[] = {"TGCAT",             // A/TGCAT, EcoT22I
                         "ATGCA"};
const char *haeIII[]  = {"CC",                // GG/CC, HaeIII
                         "GG"};
const char *hinP1I[]  = {"CGC",               // G/CGC, HinP1I
                         "GCG"};
const char *hindIII[] = {"AGCTT",             // A/AGCTT, HindIII
                         "TCGAA"};
const char *hpaII[]   = {"CGG",               // C/CGG, HpaII
                         "CCG"};
const char *kpnI[]    = {"GTACC",             // C/CATGG, KpnI
                         "GGTAC"};
const char *mluCI[]   = {"AATT",              // AATT, MluCI
                         "AATT"};
const char *mseI[]    = {"TAA",               // T/TAA, MseI
                         "TTA"};
const char *mslI[]    = {                     // CAYNN/NNRTG, MslI
        "AAATG", "AAGTG", "ACATG", "ACGTG", "AGATG", "AGGTG",
        "ATATG", "ATGTG", "CAATG", "CAGTG", "CCATG", "CCGTG",
        "CGATG", "CGGTG", "CTATG", "CTGTG", "GAATG", "GAGTG",
        "GCATG", "GCGTG", "GGATG", "GGGTG", "GTATG", "GTGTG",
        "TAATG", "TAGTG", "TCATG", "TCGTG", "TGATG", "TGGTG",
        "TTATG", "TTGTG",
        "CATTT", "CACTT", "CATGT", "CACGT", "CATCT", "CACCT",
        "CATAT", "CACAT", "CATTG", "CACTG", "CATGG", "CACGG",
        "CATCG", "CACCG", "CATAG", "CACAG", "CATTC", "CACTC",
        "CATGC", "CACGC", "CATCC", "CACCC", "CATAC", "CACAC",
        "CATTA", "CACTA", "CATGA", "CACGA", "CATCA", "CACCA",
        "CATAA", "CACAA"
};
const char *mspI[]    = {"CGG",               // C/CGG, MspI
                         "CCG"};
const char *ncoI[]    = {"CATGG",             // C/CATGG, NcoI
                         "CCATG"};
const char *ndeI[]    = {"TATG",              // CA/TATG, NdeI
                         "CATA"};
const char *nheI[]    = {"CTAGC",             // G/CTAGC, NheI
                         "GCTAG"};
const char *nlaIII[]  = {"CATG",              // CATG, NlaIII
                         "CATG"};
const char *notI[]    = {"GGCCGC",            // GC/GGCCGC, NotI
                         "GCGGCC"};
const char *nsiI[]    = {"TGCAT",             // ATGCA/T, NsiI
                         "ATGCA"};
const char *nspI[]    = {"CATGA", "CATGG",    // R/CATGY, NspI
                         "GCATG", "ACATG"};
const char *pstI[]    = {"TGCAG",             // CTGCA/G, PstI
                         "CTGCA"};
const char *rsaI[]    = {"AC",                // GT/AC, RsaI
                         "GT"};
const char *sacI[]    = {"AGCTC",             // GAGCT/C, SacI
                         "GAGCT"};
const char *sau3AI[]  = {"GATC",              // GATC, Sau3AI
                         "GATC"};
const char *sbfI[]    = {"TGCAGG",            // CCTGCA/GG, SbfI
                         "CCTGCA"};
const char *sexAI[]   = {"CCAGGT", "CCTGGT",  // A/CCWGGT, SexAI; W=A or T
                         "ACCTGG", "ACCAGG"};
const char *sgrAI[]   = {"CCGGCG", "CCGGTG",  // CR/CCGGYG, SgrAI; R=A or G; Y=C or T
                         "CGCCGG", "CACCGG"};
const char *speI[]    = {"CTAGT",             // A/CTAGT, SpeI
                         "ACTAG"};
const char *sphI[]    = {"CATGC",             // GCATG/C, SphI
                         "GCATG"};
const char *taqI[]    = {"CGA",               // T/CGA, TaqI
                         "TCG"};
const char *xbaI[]    = {"CTAGA",             // T/CTAGA, XbaI
                         "TCTAG"};
const char *xhoI[]    = {"TCGAG",             // C/TCGAG, XhoI
                         "CTCGA"};

void
initialize_renz(map<string, const char **> &renz, map<string, int> &renz_cnt, map<string, int> &renz_len) {

    renz["sbfI"]    = sbfI;    // CCTGCA/GG, SbfI
    renz["pstI"]    = pstI;    // CTGCA/G, PstI
    renz["notI"]    = notI;    // GC/GGCCGC, NotI
    renz["ecoRI"]   = ecoRI;   // G/AATTC, EcoRI
    renz["nspI"]    = nspI;
    renz["sgrAI"]   = sgrAI;   // CR/CCGGYG, SgrAI; R=A or G; Y=C or T
    renz["apeKI"]   = apeKI;   // G/CWGC, ApeKI; W=A or T
    renz["hindIII"] = hindIII; // A/AGCTT, HindIII
    renz["haeIII"]  = haeIII;  // GG/CC, HaeIII
    renz["dpnII"]   = dpnII;   // GATC, DpnII
    renz["sphI"]    = sphI;    // GCATG/C, SphI
    renz["nlaIII"]  = nlaIII;  // CATG, NlaIII
    renz["mluCI"]   = mluCI;   // AATT, MluCI
    renz["ecoT22I"] = ecoT22I; // A/TGCAT, EcoT22I
    renz["ndeI"]    = ndeI;    // CA/TATG, NdeI
    renz["nsiI"]    = nsiI;    // ATGCA/T, NsiI
    renz["mseI"]    = mseI;    // T/TAA, MseI
    renz["mslI"]    = mslI;
    renz["mspI"]    = mspI;    // C/CGG, MspI
    renz["sexAI"]   = sexAI;   // A/CCWGGT, SexAI; W=A or T
    renz["sau3AI"]  = sau3AI;  // GATC, Sau3AI
    renz["bamHI"]   = bamHI;   // G/GATCC, BamHI
    renz["xbaI"]    = xbaI;    // T/CTAGA, XbaI
    renz["eaeI"]    = eaeI;    // Y/GGCCR, EaeI
    renz["taqI"]    = taqI;    // T/CGA, TaqI
    renz["claI"]    = claI;    // AT/CGAT, ClaI
    renz["cac8I"]   = cac8I;
    renz["nheI"]    = nheI;    // G/CTAGC, NheI
    renz["speI"]    = speI;    // A/CTAGT, SpeI
    renz["apoI"]    = apoI;    // R/AATTY, ApoI, XapI
    renz["bstYI"]   = bstYI;   // R/GATCY, BstYI, PsuI
    renz["xhoI"]    = xhoI;    // C/TCGAG, XhoI
    renz["sacI"]    = sacI;    // GAGCT/C, SacI
    renz["bgIII"]   = bgIII;   // A/GATCT, BgIII
    renz["ecoRV"]   = ecoRV;   // GAT/ATC, EcoRV
    renz["kpnI"]    = kpnI;    // C/CATGG, KpnI
    renz["ddeI"]    = ddeI;    // C/TNAG, DdeI
    renz["aluI"]    = aluI;    // AG/CT, AluI
    renz["ageI"]    = ageI;    // A/CCGGT, AgeI
    renz["rsaI"]    = rsaI;    // GT/AC, RsaI
    renz["aciI"]    = aciI;    // C/CGC, AciI
    renz["bfaI"]    = bfaI;    // C/TAG, BfaI
    renz["bfuCI"]   = bfuCI;   // /GATC, BfuCI
    renz["aseI"]    = aseI;    // AT/TAAT, AseI
    renz["bspDI"]   = bspDI;   // AT/CGAT, BspDI
    renz["csp6I"]   = csp6I;   // G/TAC, Csp6I
    renz["bsaHI"]   = bsaHI;   // GR/CGYC, BsaHI
    renz["hpaII"]   = hpaII;   // C/CGG, HpaII
    renz["ncoI"]    = ncoI;    // C/CATGG, NcoI
    renz["apaLI"]   = apaLI;   // G/TGCAC, ApaLI
    renz["hinP1I"]  = hinP1I;  // G/CGC, HinP1I
    renz["bbvCI"]   = bbvCI;

    renz_cnt["sbfI"]    = 1;
    renz_cnt["pstI"]    = 1;
    renz_cnt["notI"]    = 1;
    renz_cnt["ecoRI"]   = 1;
    renz_cnt["nspI"]    = 2;
    renz_cnt["sgrAI"]   = 2;
    renz_cnt["apeKI"]   = 2;
    renz_cnt["hindIII"] = 1;
    renz_cnt["haeIII"]  = 1;
    renz_cnt["dpnII"]   = 1;
    renz_cnt["sphI"]    = 1;
    renz_cnt["nlaIII"]  = 1;
    renz_cnt["mluCI"]   = 1;
    renz_cnt["ecoT22I"] = 1;
    renz_cnt["ndeI"]    = 1;
    renz_cnt["nsiI"]    = 1;
    renz_cnt["mseI"]    = 1;
    renz_cnt["mslI"]    = 32;
    renz_cnt["mspI"]    = 1;
    renz_cnt["sexAI"]   = 2;
    renz_cnt["sau3AI"]  = 1;
    renz_cnt["bamHI"]   = 1;
    renz_cnt["xbaI"]    = 1;
    renz_cnt["eaeI"]    = 2;
    renz_cnt["taqI"]    = 1;
    renz_cnt["claI"]    = 1;
    renz_cnt["cac8I"]   = 4;
    renz_cnt["nheI"]    = 1;
    renz_cnt["speI"]    = 1;
    renz_cnt["apoI"]    = 2;
    renz_cnt["bstYI"]   = 2;
    renz_cnt["xhoI"]    = 1;
    renz_cnt["sacI"]    = 1;
    renz_cnt["bgIII"]   = 1;
    renz_cnt["ecoRV"]   = 1;
    renz_cnt["kpnI"]    = 1;
    renz_cnt["ddeI"]    = 4;
    renz_cnt["aluI"]    = 1;
    renz_cnt["ageI"]    = 1;
    renz_cnt["rsaI"]    = 1;
    renz_cnt["aciI"]    = 2;
    renz_cnt["bfaI"]    = 1;
    renz_cnt["bfuCI"]   = 1;
    renz_cnt["aseI"]    = 1;
    renz_cnt["bspDI"]   = 1;
    renz_cnt["csp6I"]   = 1;
    renz_cnt["bsaHI"]   = 2;
    renz_cnt["hpaII"]   = 1;
    renz_cnt["ncoI"]    = 1;
    renz_cnt["apaLI"]   = 1;
    renz_cnt["hinP1I"]  = 1;
    renz_cnt["bbvCI"]   = 1;

    renz_len["sbfI"]    = 6;
    renz_len["pstI"]    = 5;
    renz_len["notI"]    = 6;
    renz_len["ecoRI"]   = 5;
    renz_len["nspI"]    = 5;
    renz_len["sgrAI"]   = 6;
    renz_len["apeKI"]   = 4;
    renz_len["hindIII"] = 5;
    renz_len["haeIII"]  = 2;
    renz_len["dpnII"]   = 4;
    renz_len["sphI"]    = 5;
    renz_len["nlaIII"]  = 4;
    renz_len["mluCI"]   = 4;
    renz_len["ecoT22I"] = 5;
    renz_len["ndeI"]    = 4;
    renz_len["nsiI"]    = 5;
    renz_len["mseI"]    = 3;
    renz_len["mslI"]    = 5;
    renz_len["mspI"]    = 3;
    renz_len["sexAI"]   = 6;
    renz_len["sau3AI"]  = 4;
    renz_len["bamHI"]   = 5;
    renz_len["xbaI"]    = 5;
    renz_len["eaeI"]    = 5;
    renz_len["taqI"]    = 3;
    renz_len["claI"]    = 4;
    renz_len["cac8I"]   = 3;
    renz_len["nheI"]    = 5;
    renz_len["speI"]    = 5;
    renz_len["apoI"]    = 5;
    renz_len["bstYI"]   = 5;
    renz_len["xhoI"]    = 5;
    renz_len["sacI"]    = 5;
    renz_len["bgIII"]   = 5;
    renz_len["ecoRV"]   = 3;
    renz_len["kpnI"]    = 5;
    renz_len["ddeI"]    = 4;
    renz_len["aluI"]    = 2;
    renz_len["ageI"]    = 5;
    renz_len["rsaI"]    = 2;
    renz_len["aciI"]    = 3;
    renz_len["bfaI"]    = 3;
    renz_len["bfuCI"]   = 4;
    renz_len["aseI"]    = 4;
    renz_len["bspDI"]   = 4;
    renz_len["csp6I"]   = 3;
    renz_len["bsaHI"]   = 4;
    renz_len["hpaII"]   = 3;
    renz_len["ncoI"]    = 5;
    renz_len["apaLI"]   = 5;
    renz_len["hinP1I"]  = 3;
    renz_len["bbvCI"]   = 5;
}

void
initialize_renz_olap(map<string, int> &renz_olap) {
    renz_olap["sbfI"]   = 4;
}

