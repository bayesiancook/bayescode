#ifndef BIOLOGICALSEQUENCES_H
#define BIOLOGICALSEQUENCES_H

#include <cstdlib>
#include <iostream>
#include <string>

enum DataType { DNA = 0, RNA = 1, Protein = 2, Other = 3 };

const int RYN = 2;
const char RYset[] = {'R', 'Y'};

const char AAset[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                      'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'a', 'c', 'd', 'e', 'f', 'g',
                      'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w',
                      'y', '-', '?', '$', '.', 'B', 'Z', '*', 'X', 'x'};
const int AAN = 49;
const int DNAN = 37;
const int RNAN = 37;
const char DNAset[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'B', 'D', 'H', 'K', 'M',
                       'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                       's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};
const char RNAset[] = {'A', 'C', 'G', 'U', 'a', 'c', 'g', 'u', 'B', 'D', 'H', 'K', 'M',
                       'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                       's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};

// amino acids

const int precision = 10000;
const std::string Alphabet = "Amino_Acids";
const char AminoAcids[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                           'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'};
const char aminoacids[] = {'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
                           'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', '-'};
const char RYletters[] = {'R', 'Y'};
const char DNAletters[] = {'A', 'C', 'G', 'T'};
const char dnaletters[] = {'a', 'c', 'g', 't'};
const char RNAletters[] = {'A', 'C', 'G', 'U'};
const char rnaletters[] = {'a', 'c', 'g', 'u'};

const int Dayhoff6Table[] = {3, 5, 2, 2, 4, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 4, 4};
const int Dayhoff4Table[] = {3, -1, 2, 2, 0, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 2, 2};

const int unknown = -1;

const int Naa = 20;
const int Naarr = Naa * (Naa-1) / 2;
const int Nnuc = 4;
const int Nrr = Nnuc * (Nnuc - 1) / 2;
const int Ncodon = 64;
const std::string Codons[] = {
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT",
    "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
    "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA",
    "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG",
    "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"};
const int codonpos[][64] = {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                            {3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1,
                             1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0,
                             2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2},
                            {3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1,
                             0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2,
                             3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2}};

enum GeneticCodeType { Universal = 0, MtMam = 1, MtInv = 2, MtProt = 3, MtEch = 4 };
// universal genetic code
// const std::string UniStopCodons[] = {"TAA","TAG","TGA"};
const int UniNStopCodons = 3;
const int UniStopCodons[] = {10, 11, 14};
int const UniCodonCode[] = {4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  -1, 18,
                            9,  9,  9,  9,  12, 12, 12, 12, 6,  6,  13, 13, 14, 14, 14, 14,
                            7,  7,  7,  10, 16, 16, 16, 16, 11, 11, 8,  8,  15, 15, 14, 14,
                            17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const int UniStopPos1[] = {3, 3, 3};
const int UniStopPos2[] = {0, 0, 2};
const int UniStopPos3[] = {0, 2, 0};

const int MtInvNStopCodons = 2;
const int MtInvStopCodons[] = {10, 11};
int const MtInvCodonCode[] = {4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18,
                              9,  9,  9,  9,  12, 12, 12, 12, 6,  6,  13, 13, 14, 14, 14, 14,
                              7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,  15, 15, 15, 15,
                              17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const int MtInvStopPos1[] = {3, 3};
const int MtInvStopPos2[] = {0, 0};
const int MtInvStopPos3[] = {0, 2};

// mammal mitochondrial genetic code
const int MtMamNStopCodons = 4;
const int MtMamStopCodons[] = {10, 11, 46, 47};
int const MtMamCodonCode[] = {4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18,
                              9,  9,  9,  9,  12, 12, 12, 12, 6,  6,  13, 13, 14, 14, 14, 14,
                              7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,  15, 15, -1, -1,
                              17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const int MtMamStopPos1[] = {3, 3, 0, 0};
const int MtMamStopPos2[] = {0, 0, 2, 2};
const int MtMamStopPos3[] = {0, 2, 0, 2};

/*
// Protozoan and Coelenterate mitochondrial genetic code
const int MtProtNStopCodons = 2;
const int MtProtStopCodons[] = {10,11};
int const MtProtCodonCode[] =
{4,4,9,9,15,15,15,15,19,19,-1,-1,1,1,18,18,9,9,9,9,12,12,12,12,6,6,13,13,14,14,14,14,7,7,7,10,16,16,16,16,11,11,8,8,15,15,14,14,17,17,17,17,0,0,0,0,2,2,3,3,5,5,5,5};

// Echinoderm and Flatworm mitochondrial genetic code
const int MtEchNStopCodons = 2;
const int MtEchStopCodons[] = {10,11};
int const MtEchCodonCode[] =
{4,4,9,9,15,15,15,15,19,19,-1,-1,1,1,18,18,9,9,9,9,12,12,12,12,6,6,13,13,14,14,14,14,7,7,7,10,16,16,16,16,11,11,11,8,15,15,15,15,17,17,17,17,0,0,0,0,2,2,3,3,5,5,5,5};
*/

inline std::istream &operator>>(std::istream &is, GeneticCodeType &type) {
    std::string t;
    is >> t;
    if (t == "Universal") {
        type = Universal;
    } else if (t == "MtMam") {
        type = MtMam;
    } else if (t == "MtInv") {
        type = MtInv;
    } else if (t == "MtProt") {
        type = MtProt;
    } else if (t == "MtEch") {
        type = MtEch;
    } else {
        std::cerr << "error in std::istream genetic code type\n";
        std::cerr << type << '\n';
        exit(1);
    }
    return is;
}

inline std::ostream &operator<<(std::ostream &os, GeneticCodeType type) {
    if (type == Universal) {
        os << "Universal\n";
    } else if (type == MtMam) {
        os << "MtMam\n";
    } else if (type == MtInv) {
        os << "MtInv\n";
    } else if (type == MtProt) {
        os << "MtProt\n";
    } else if (type == MtEch) {
        os << "MtEch\n";
    } else {
        std::cerr << "error in std::ostream genetic code type\n";
        std::cerr << (int)type << '\n';
        exit(1);
    }
    return os;
}

const double LG_RR[] = {2.43501,0.386559,1.01598,0.248188,2.02114,0.35106,0.146574,0.524859,0.386747,1.09961,0.270803,1.15206,0.94882,0.415855,4.62446,2.09301,2.49251,0.176789,0.214201,0.0611966,0.00342322,1.08123,0.556892,0.626623,0.313658,0.0129776,0.581098,0.874258,0.517278,0.0737435,0.0829654,0.522935,2.72396,1.11863,1.91672,0.655562,1.1402,5.12993,0.0170374,0.826565,0.906972,0.010458,0.27681,0.014748,0.0249932,4.96585,0.385885,0.512014,0.121261,1.21333,0.416606,0.0371423,0.0292405,0.132172,0.0184023,0.341267,0.414673,0.0433031,1.7679,0.0681587,0.16996,0.529941,0.410296,4.03888,0.356061,0.598676,0.591407,0.23971,0.0761602,0.117429,0.0876389,0.667316,1.08855,0.0233983,2.53635,1.75976,0.0875798,0.0924114,0.035076,0.0515759,0.353957,0.161416,0.640457,2.40372,7.63435,0.304716,0.00851629,0.290191,0.0432987,0.136505,1.40641,0.192681,0.262136,0.381714,1.70218,0.127015,0.0750342,0.262657,0.0534908,0.106516,0.682111,0.358356,0.432853,4.41125,0.497794,4.70891,2.37387,0.968501,0.571566,0.116427,0.584078,5.19151,0.155611,4.055,4.18072,0.187342,0.0765801,0.0712711,0.124231,0.0627125,1.01127,10.4177,0.109234,0.22747,0.134512,0.642335,2.09847,0.381839,3.16401,6.1886,0.732413,1.11216,0.181177,0.048821,0.129065,6.17517,0.0669398,0.243648,0.5698,0.295289,0.178326,0.296353,1.66574,0.606165,0.293136,0.362942,0.0976794,1.63622,0.473613,0.339422,1.97647,1.85746,0.681046,0.470848,0.158271,1.6589,0.73554,3.92125,1.95721,0.081869,0.0443899,0.598725,0.610728,0.325308,1.30906,0.559051,0.290058,0.0930628,0.0876663,2.74689,1.19724,1.05666,0.205761,0.231066,0.251745,0.839504,0.566405,0.167174,0.580707,0.307606,6.33164,0.096231,0.243453,0.391844,2.14061,0.137765,0.240498,0.185392,0.243896,3.08335};

const double LG_Stat[] = {
0.079066,
0.0129369,
0.0530516,
0.0715863,
0.0423017,
0.0573372,
0.0223546,
0.0621565,
0.0646003,
0.099081,
0.0229506,
0.0419774,
0.0440395,
0.0407668,
0.0559413,
0.0611971,
0.0532871,
0.0691469,
0.0120656,
0.0341554
};

#endif  // BIOLOGICALSEQUENCES_H
