#include "genotype_dictionaries.h"

int
encoded_gtypes[4][4] =
{
     // A  C  G   T
    {1, 2, 3,  4}, // A
    {2, 5, 6,  7}, // C
    {3, 6, 8,  9}, // G
    {4, 7, 9, 10}  // T
};

void
initialize_dictionaries(map<string, map<string, string>> &global_dictionary)
{
    global_dictionary["ab/--"]["aa"] = "aa";
    global_dictionary["ab/--"]["bb"] = "bb";

    global_dictionary["--/ab"]["aa"] = "aa";
    global_dictionary["--/ab"]["bb"] = "bb";

    global_dictionary["aa/bb"]["aa"] = "aa";
    global_dictionary["aa/bb"]["ab"] = "ab";
    global_dictionary["aa/bb"]["bb"] = "bb";

    global_dictionary["ab/ac"]["aa"] = "aa";
    global_dictionary["ab/ac"]["ab"] = "ab";
    global_dictionary["ab/ac"]["bb"] = "bb";
    global_dictionary["ab/ac"]["ac"] = "ac";
    global_dictionary["ab/ac"]["cc"] = "cc";
    global_dictionary["ab/ac"]["bc"] = "bc";

    global_dictionary["ab/cd"]["aa"] = "aa";
    global_dictionary["ab/cd"]["ab"] = "ab";
    global_dictionary["ab/cd"]["bb"] = "bb";
    global_dictionary["ab/cd"]["cc"] = "cc";
    global_dictionary["ab/cd"]["cd"] = "cd";
    global_dictionary["ab/cd"]["dd"] = "dd";
    global_dictionary["ab/cd"]["ac"] = "ac";
    global_dictionary["ab/cd"]["ad"] = "ad";
    global_dictionary["ab/cd"]["bc"] = "bc";
    global_dictionary["ab/cd"]["bd"] = "bd";

    global_dictionary["ab/aa"]["aa"] = "aa";
    global_dictionary["ab/aa"]["ab"] = "ab";
    global_dictionary["ab/aa"]["bb"] = "bb";

    global_dictionary["aa/ab"]["aa"] = "aa";
    global_dictionary["aa/ab"]["ab"] = "ab";
    global_dictionary["aa/ab"]["bb"] = "bb";

    global_dictionary["ab/cc"]["aa"] = "aa";
    global_dictionary["ab/cc"]["ab"] = "ab";
    global_dictionary["ab/cc"]["bb"] = "bb";
    global_dictionary["ab/cc"]["cc"] = "cc";
    global_dictionary["ab/cc"]["ac"] = "ac";
    global_dictionary["ab/cc"]["bc"] = "bc";

    global_dictionary["cc/ab"]["aa"] = "aa";
    global_dictionary["cc/ab"]["ab"] = "ab";
    global_dictionary["cc/ab"]["bb"] = "bb";
    global_dictionary["cc/ab"]["cc"] = "cc";
    global_dictionary["cc/ab"]["ac"] = "ac";
    global_dictionary["cc/ab"]["bc"] = "bc";

    global_dictionary["ab/ab"]["aa"] = "aa";
    global_dictionary["ab/ab"]["bb"] = "bb";
    global_dictionary["ab/ab"]["ab"] = "ab";
}

void
load_cp_dictionary(map<string, string> &types,
                   map<string, map<string, string>> &dictionary)
{
    types["ab/--"] = "ab/--";
    types["--/ab"] = "--/ab";
    types["ab/aa"] = "ab/aa";
    types["aa/ab"] = "aa/ab";
    types["ab/a-"] = "ab/--";
    types["-a/ab"] = "--/ab";
    types["ab/c-"] = "ab/cd";
    types["-c/ab"] = "ab/cd";
    types["ab/cc"] = "ab/--";
    types["cc/ab"] = "--/ab";
    types["ab/ab"] = "ab/ab";
    types["ab/ac"] = "ab/ac";
    types["ab/cd"] = "ab/cd";
    types["-a/bb"] = "ab/--";
    types["aa/b-"] = "--/ab";

    dictionary["ab/--"]["--"] = "--";
    dictionary["ab/--"]["aa"] = "aa";
    dictionary["ab/--"]["ab"] = "ab";
    dictionary["ab/--"]["bb"] = "ab";
    dictionary["ab/--"]["ac"] = "aa";
    dictionary["ab/--"]["bc"] = "ab";

    dictionary["--/ab"]["--"] = "--";
    dictionary["--/ab"]["aa"] = "aa";
    dictionary["--/ab"]["ab"] = "ab";
    dictionary["--/ab"]["bb"] = "ab";
    dictionary["--/ab"]["ac"] = "aa";
    dictionary["--/ab"]["bc"] = "ab";

    dictionary["ab/aa"]["--"] = "--";
    dictionary["ab/aa"]["aa"] = "aa";
    dictionary["ab/aa"]["ab"] = "ab";

    dictionary["aa/ab"]["--"] = "--";
    dictionary["aa/ab"]["aa"] = "aa";
    dictionary["aa/ab"]["ab"] = "ab";

    dictionary["ab/ab"]["--"] = "--";
    dictionary["ab/ab"]["ab"] = "ab";
    dictionary["ab/ab"]["aa"] = "aa";
    dictionary["ab/ab"]["bb"] = "bb";

    dictionary["ab/ac"]["--"] = "--";
    dictionary["ab/ac"]["ab"] = "ab";
    dictionary["ab/ac"]["ac"] = "ac";
    dictionary["ab/ac"]["bc"] = "bc";
    dictionary["ab/ac"]["aa"] = "aa";

    dictionary["ab/cd"]["--"] = "--";
    dictionary["ab/cd"]["ac"] = "ac";
    dictionary["ab/cd"]["ad"] = "ad";
    dictionary["ab/cd"]["bc"] = "bc";
    dictionary["ab/cd"]["bd"] = "bd";

    return;
}

void
load_joinmap_cp_dictionary(map<string, string> &types,
                           map<string, map<string, string>> &dictionary)
{
    types["ab/--"] = "lmx--";
    types["--/ab"] = "--xnp";
    types["ab/aa"] = "lmxll";
    types["aa/ab"] = "nnxnp";
    types["ab/ab"] = "hkxhk";
    types["ab/ac"] = "efxeg";
    types["ab/cd"] = "abxcd";

    dictionary["lmx--"]["--"] = "--";
    dictionary["lmx--"]["aa"] = "ll";
    dictionary["lmx--"]["ab"] = "lm";
    dictionary["lmx--"]["bb"] = "lm";
    dictionary["lmx--"]["ac"] = "ll";
    dictionary["lmx--"]["bc"] = "lm";

    dictionary["--xnp"]["--"] = "--";
    dictionary["--xnp"]["aa"] = "nn";
    dictionary["--xnp"]["ab"] = "np";
    dictionary["--xnp"]["bb"] = "np";
    dictionary["--xnp"]["ac"] = "nn";
    dictionary["--xnp"]["bc"] = "np";

    dictionary["lmxll"]["--"] = "--";
    dictionary["lmxll"]["aa"] = "ll";
    dictionary["lmxll"]["ab"] = "lm";

    dictionary["nnxnp"]["--"] = "--";
    dictionary["nnxnp"]["aa"] = "nn";
    dictionary["nnxnp"]["ab"] = "np";

    dictionary["hkxhk"]["--"] = "--";
    dictionary["hkxhk"]["ab"] = "hk";
    dictionary["hkxhk"]["aa"] = "hh";
    dictionary["hkxhk"]["bb"] = "kk";

    dictionary["efxeg"]["--"] = "--";
    dictionary["efxeg"]["ab"] = "ef";
    dictionary["efxeg"]["ac"] = "eg";
    dictionary["efxeg"]["bc"] = "fg";
    dictionary["efxeg"]["aa"] = "ee";

    dictionary["abxcd"]["--"] = "--";
    dictionary["abxcd"]["ac"] = "ac";
    dictionary["abxcd"]["ad"] = "ad";
    dictionary["abxcd"]["bc"] = "bc";
    dictionary["abxcd"]["bd"] = "bd";

    return;
}

void
load_onemap_cp_dictionary(map<string, string> &types,
                          map<string, map<string, string>> &dictionary)
{
    types["ab/--"] = "abxoo";
    types["--/ab"] = "ooxab";
    types["ab/aa"] = "abxaa";
    types["aa/ab"] = "aaxab";
    types["ab/ab"] = "abxab";
    types["ab/ac"] = "abxac";
    types["ab/cd"] = "abxcd";

    // D1.11
    dictionary["abxoo"]["--"] = "-";
    dictionary["abxoo"]["aa"] = "a";
    dictionary["abxoo"]["bb"] = "b";

    // D2.16
    dictionary["ooxab"]["--"] = "-";
    dictionary["ooxab"]["aa"] = "a";
    dictionary["ooxab"]["bb"] = "b";

    // D1.10
    dictionary["abxaa"]["--"] = "-";
    dictionary["abxaa"]["aa"] = "a";
    dictionary["abxaa"]["ab"] = "ab";

    // D2.15
    dictionary["aaxab"]["--"] = "-";
    dictionary["aaxab"]["aa"] = "a";
    dictionary["aaxab"]["ab"] = "ab";

    // B3.7
    dictionary["abxab"]["--"] = "-";
    dictionary["abxab"]["ab"] = "ab";
    dictionary["abxab"]["aa"] = "a";
    dictionary["abxab"]["bb"] = "b";

    // A.2
    dictionary["abxac"]["--"] = "-";
    dictionary["abxac"]["ab"] = "ba";
    dictionary["abxac"]["ac"] = "ac";
    dictionary["abxac"]["bc"] = "bc";
    dictionary["abxac"]["aa"] = "a";

    // A.1
    dictionary["abxcd"]["--"] = "-";
    dictionary["abxcd"]["ac"] = "ac";
    dictionary["abxcd"]["ad"] = "ad";
    dictionary["abxcd"]["bc"] = "bc";
    dictionary["abxcd"]["bd"] = "bd";

    return;
}

void
load_bc_dictionary(map<string, string> &types,
                   map<string, map<string, string>>  &dictionary)
{
    types["aa/bb"] = "aa/bb";
    types["bb/aa"] = "bb/aa";
    types["ab/cc"] = "ab/cc";
    types["cc/ab"] = "cc/ab";

    dictionary["aa/bb"]["--"] = "--";
    dictionary["aa/bb"]["aa"] = "aa";
    dictionary["aa/bb"]["ab"] = "ab";
    dictionary["aa/bb"]["bb"] = "bb";

    dictionary["bb/aa"]["--"] = "--";
    dictionary["bb/aa"]["aa"] = "aa";
    dictionary["bb/aa"]["ab"] = "ab";
    dictionary["bb/aa"]["bb"] = "bb";

    dictionary["ab/cc"]["--"] = "--";
    dictionary["ab/cc"]["ac"] = "ac";
    dictionary["ab/cc"]["bc"] = "bc";
    dictionary["ab/cc"]["ab"] = "ab";
    dictionary["ab/cc"]["aa"] = "aa";
    dictionary["ab/cc"]["bb"] = "bb";

    dictionary["cc/ab"]["--"] = "--";
    dictionary["cc/ab"]["ac"] = "ac";
    dictionary["cc/ab"]["bc"] = "bc";
    dictionary["cc/ab"]["ab"] = "ab";
    dictionary["cc/ab"]["aa"] = "aa";
    dictionary["cc/ab"]["bb"] = "bb";
}

void
load_f2_dictionary(map<string, string> &types,
                   map<string, map<string, string>> &dictionary)
{
    types["aa/bb"] = "aa/bb";
    types["ab/cd"] = "ab/cd";
    types["ab/aa"] = "ab/aa";
    types["aa/ab"] = "aa/ab";
    types["ab/cc"] = "ab/cc";
    types["cc/ab"] = "cc/ab";

    dictionary["aa/bb"]["aa"] = "aa";
    dictionary["aa/bb"]["ab"] = "ab";
    dictionary["aa/bb"]["bb"] = "bb";
    dictionary["aa/bb"]["--"] = "--";

    dictionary["ab/cd"]["aa"] = "aa";
    dictionary["ab/cd"]["ab"] = "ab";
    dictionary["ab/cd"]["bb"] = "bb";
    dictionary["ab/cd"]["cc"] = "cc";
    dictionary["ab/cd"]["cd"] = "cd";
    dictionary["ab/cd"]["dd"] = "dd";
    dictionary["ab/cd"]["ac"] = "ac";
    dictionary["ab/cd"]["ad"] = "ad";
    dictionary["ab/cd"]["bc"] = "bc";
    dictionary["ab/cd"]["bd"] = "bd";
    dictionary["ab/cd"]["--"] = "--";

    dictionary["ab/aa"]["aa"] = "--";
    dictionary["ab/aa"]["ab"] = "--";
    dictionary["ab/aa"]["bb"] = "bb";
    dictionary["ab/aa"]["--"] = "--";

    dictionary["aa/ab"]["aa"] = "--";
    dictionary["aa/ab"]["ab"] = "--";
    dictionary["aa/ab"]["bb"] = "bb";
    dictionary["aa/ab"]["--"] = "--";

    dictionary["ab/cc"]["aa"] = "aa";
    dictionary["ab/cc"]["ab"] = "ab";
    dictionary["ab/cc"]["bb"] = "bb";
    dictionary["ab/cc"]["cc"] = "cc";
    dictionary["ab/cc"]["ac"] = "--";
    dictionary["ab/cc"]["bc"] = "--";
    dictionary["ab/cc"]["--"] = "--";

    dictionary["cc/ab"]["aa"] = "aa";
    dictionary["cc/ab"]["ab"] = "ab";
    dictionary["cc/ab"]["bb"] = "bb";
    dictionary["cc/ab"]["cc"] = "cc";
    dictionary["cc/ab"]["ac"] = "--";
    dictionary["cc/ab"]["bc"] = "--";
    dictionary["cc/ab"]["--"] = "--";
}

void
load_mm_bc_dictionary(map<string, string> &types,
                      map<string, map<string, string>> &dictionary)
{
    types["aa/bb"] = "aaxbb";
    types["bb/aa"] = "bbxaa";
    types["ab/cc"] = "abxcc";
    types["cc/ab"] = "ccxab";

    dictionary["aaxbb"]["--"] = "-";
    dictionary["aaxbb"]["aa"] = "b";
    dictionary["aaxbb"]["ab"] = "h";
    dictionary["aaxbb"]["bb"] = "h";

    dictionary["bbxaa"]["--"] = "-";
    dictionary["bbxaa"]["aa"] = "h";
    dictionary["bbxaa"]["ab"] = "h";
    dictionary["bbxaa"]["bb"] = "a";

    dictionary["abxcc"]["--"] = "-";
    dictionary["abxcc"]["ac"] = "h";
    dictionary["abxcc"]["bc"] = "h";
    dictionary["abxcc"]["ab"] = "b";
    dictionary["abxcc"]["aa"] = "b";
    dictionary["abxcc"]["bb"] = "b";

    dictionary["ccxab"]["--"] = "-";
    dictionary["ccxab"]["ac"] = "h";
    dictionary["ccxab"]["bc"] = "h";
    dictionary["ccxab"]["ab"] = "a";
    dictionary["ccxab"]["aa"] = "a";
    dictionary["ccxab"]["bb"] = "a";
}

void
load_mm_f2_dictionary(map<string, string> &types,
                      map<string, map<string, string>> &dictionary)
{
    types["aa/bb"] = "aaxbb";
    types["ab/cd"] = "abxcd";
    types["ab/aa"] = "abxaa";
    types["aa/ab"] = "aaxab";
    types["ab/cc"] = "abxcc";
    types["cc/ab"] = "ccxab";

    dictionary["aaxbb"]["aa"] = "a";
    dictionary["aaxbb"]["ab"] = "h";
    dictionary["aaxbb"]["bb"] = "b";
    dictionary["aaxbb"]["--"] = "-";

    dictionary["abxcd"]["aa"] = "a";
    dictionary["abxcd"]["ab"] = "a";
    dictionary["abxcd"]["bb"] = "a";
    dictionary["abxcd"]["cc"] = "b";
    dictionary["abxcd"]["cd"] = "b";
    dictionary["abxcd"]["dd"] = "b";
    dictionary["abxcd"]["ac"] = "h";
    dictionary["abxcd"]["ad"] = "h";
    dictionary["abxcd"]["bc"] = "h";
    dictionary["abxcd"]["bd"] = "h";
    dictionary["abxcd"]["--"] = "-";

    dictionary["abxaa"]["aa"] = "-";
    dictionary["abxaa"]["ab"] = "-";
    dictionary["abxaa"]["bb"] = "a";
    dictionary["abxaa"]["--"] = "-";

    dictionary["aaxab"]["aa"] = "-";
    dictionary["aaxab"]["ab"] = "-";
    dictionary["aaxab"]["bb"] = "b";
    dictionary["aaxab"]["--"] = "-";

    dictionary["abxcc"]["aa"] = "a";
    dictionary["abxcc"]["ab"] = "a";
    dictionary["abxcc"]["bb"] = "a";
    dictionary["abxcc"]["cc"] = "b";
    dictionary["abxcc"]["ac"] = "-";
    dictionary["abxcc"]["bc"] = "-";
    dictionary["abxcc"]["--"] = "-";

    dictionary["ccxab"]["aa"] = "b";
    dictionary["ccxab"]["ab"] = "b";
    dictionary["ccxab"]["bb"] = "b";
    dictionary["ccxab"]["cc"] = "a";
    dictionary["ccxab"]["ac"] = "-";
    dictionary["ccxab"]["bc"] = "-";
    dictionary["ccxab"]["--"] = "-";
}

void
load_dh_dictionary(map<string, string> &types,
                   map<string, map<string, string>> &dictionary)
{
    types["ab/--"] = "ab/--";
    types["--/ab"] = "--/ab";

    dictionary["ab/--"]["aa"] = "aa";
    dictionary["ab/--"]["bb"] = "bb";
    dictionary["ab/--"]["--"] = "--";

    dictionary["--/ab"]["aa"] = "aa";
    dictionary["--/ab"]["bb"] = "bb";
    dictionary["--/ab"]["--"] = "--";
}

void
load_mm_dh_dictionary(map<string, string> &types,
                      map<string, map<string, string>> &dictionary)
{
    types["ab/--"] = "abx--";
    types["--/ab"] = "--xab";

    dictionary["abx--"]["aa"] = "a";
    dictionary["abx--"]["bb"] = "b";
    dictionary["abx--"]["--"] = "-";

    dictionary["--xab"]["aa"] = "a";
    dictionary["--xab"]["bb"] = "b";
    dictionary["--xab"]["--"] = "-";
}

void
load_segregation_ratios(CrossT type,
                        map<string, map<string, double>> &segregation_ratios)
{
    switch(type) {
    case CrossT::cp:
        segregation_ratios["ab/--"]["aa"] = 0.50;
        segregation_ratios["ab/--"]["ab"] = 0.50;

        segregation_ratios["--/ab"]["aa"] = 0.50;
        segregation_ratios["--/ab"]["ab"] = 0.50;

        segregation_ratios["ab/aa"]["aa"] = 0.50;
        segregation_ratios["ab/aa"]["ab"] = 0.50;

        segregation_ratios["aa/ab"]["aa"] = 0.50;
        segregation_ratios["aa/ab"]["ab"] = 0.50;

        segregation_ratios["ab/ab"]["ab"] = 0.50;
        segregation_ratios["ab/ab"]["aa"] = 0.25;
        segregation_ratios["ab/ab"]["bb"] = 0.25;

        segregation_ratios["ab/ac"]["ab"] = 0.25;
        segregation_ratios["ab/ac"]["ac"] = 0.25;
        segregation_ratios["ab/ac"]["bc"] = 0.25;
        segregation_ratios["ab/ac"]["aa"] = 0.25;

        segregation_ratios["ab/cd"]["ac"] = 0.25;
        segregation_ratios["ab/cd"]["ad"] = 0.25;
        segregation_ratios["ab/cd"]["bc"] = 0.25;
        segregation_ratios["ab/cd"]["bd"] = 0.25;
        break;
    case CrossT::f2:
        segregation_ratios["aa/bb"]["a"] = 0.25;
        segregation_ratios["aa/bb"]["b"] = 0.25;
        segregation_ratios["aa/bb"]["h"] = 0.50;

        segregation_ratios["ab/cd"]["a"] = 0.25;
        segregation_ratios["ab/cd"]["b"] = 0.25;
        segregation_ratios["ab/cd"]["h"] = 0.50;

        segregation_ratios["ab/aa"]["a"] = 1.00;

        segregation_ratios["aa/ab"]["b"] = 1.00;

        segregation_ratios["ab/cc"]["a"] = 0.50;
        segregation_ratios["ab/cc"]["b"] = 0.50;

        segregation_ratios["cc/ab"]["b"] = 0.50;
        segregation_ratios["cc/ab"]["a"] = 0.50;
        break;
    case CrossT::bc1:
        segregation_ratios["aa/bb"]["h"] = 0.50;
        segregation_ratios["aa/bb"]["b"] = 0.50;

        segregation_ratios["bb/aa"]["h"] = 0.50;
        segregation_ratios["bb/aa"]["a"] = 0.50;

        segregation_ratios["ab/cc"]["h"] = 0.50;
        segregation_ratios["ab/cc"]["b"] = 0.50;

        segregation_ratios["cc/ab"]["h"] = 0.50;
        segregation_ratios["cc/ab"]["a"] = 0.50;
        break;
    case CrossT::dh:
        segregation_ratios["ab/--"]["a"] = 0.50;
        segregation_ratios["ab/--"]["b"] = 0.50;

        segregation_ratios["--/ab"]["a"] = 0.50;
        segregation_ratios["--/ab"]["b"] = 0.50;
    break;
    case CrossT::unk:
    default:
        break;
    }

    return;
}
