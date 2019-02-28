#include "CodonSequenceAlignment.hpp"
#include <cstdlib>
#include <iostream>
#include "Exception.hpp"
#include "Random.hpp"
using namespace std;

CodonSequenceAlignment::CodonSequenceAlignment(
    SequenceAlignment *from, bool force_stops, GeneticCodeType type) {
    try {
        DNAsource = from;

        if (from->GetNsite() % 3 != 0) {
            cerr << "not multiple of three\n";
            exit(1);
        }
        Nsite = from->GetNsite() / 3;
        Ntaxa = from->GetNtaxa();
        auto tempstatespace = new CodonStateSpace(type);
        statespace = tempstatespace;

        taxset = DNAsource->GetTaxonSet();
        owntaxset = false;

        // make my own arrays
        // make translation
        Data.assign(Ntaxa, std::vector<int>(Nsite, 0));
        for (int i = 0; i < Ntaxa; i++) {
            for (int j = 0; j < Nsite; j++) {
                try {
                    Data[i][j] =
                        GetCodonStateSpace()->GetCodonFromDNA(DNAsource->GetState(i, 3 * j),
                            DNAsource->GetState(i, 3 * j + 1), DNAsource->GetState(i, 3 * j + 2));
                    if (Data[i][j] == -1) {
                        if ((DNAsource->GetState(i, 3 * j) != -1) &&
                            (DNAsource->GetState(i, 3 * j + 1) != -1) &&
                            (DNAsource->GetState(i, 3 * j + 2) != -1)) {
                            // cerr << "in CodonSequenceAlignment: taxon " <<
                            // taxset->GetTaxon(i) <<
                            // " and codon " << j+1 << " (site " << 3*j+1 << ") :";
                            // cerr << nucspace->GetState(DNAsource->GetState(i, 3*j)) <<
                            // nucspace->GetState(DNAsource->GetState(i, 3*j+1)) <<
                            // nucspace->GetState(DNAsource->GetState(i, 3*j+2)) << '\n';
                        }
                    }
                } catch (...) {
                    // catch(Exception e)	{
                    // cerr << "in CodonSequenceAlignment: taxon " << i << " and codon "
                    // << j << "
                    // (site " << 3*j << ")\n";
                    // cerr << "taxon : " << taxset->GetTaxon(i) << '\n';
                    if (force_stops) {
                        // Data[i][j] = -2;
                        Data[i][j] = -1;
                    } else {
                        throw;
                    }
                }
            }
        }

    } catch (Exception) {
        cerr << "Codon Sequence Alignment: failed to read the datafile\n";
        exit(1);
    }
}

void CodonSequenceAlignment::ToStream(ostream &os) {
    os << Ntaxa << '\t' << 3 * Nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) { max = l; }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) { os << ' '; }
        for (int j = 0; j < Nsite; j++) { os << statespace->GetState(GetState(i, j)); }
        os << '\n';
    }
    os << '\n';
}

void CodonSequenceAlignment::ToStream(ostream &os, int pos) {
    int included = 0;
    for (int i = 0; i < Nsite; i++) {
        if (!AllMissingColumn(i)) { included++; }
    }
    os << Ntaxa << '\t' << included << '\n';
    // os << Ntaxa << '\t' << Nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) { max = l; }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) { os << ' '; }
        for (int j = 0; j < Nsite; j++) {
            if (!AllMissingColumn(j)) {
                os << GetCodonStateSpace()->GetDNAStateSpace()->GetState(
                    GetCodonStateSpace()->GetCodonPosition(pos, GetState(i, j)));
            }
        }
        os << '\n';
    }
    os << '\n';
}
