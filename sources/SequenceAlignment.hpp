#ifndef SEQUENCEALIGNMENT_H
#define SEQUENCEALIGNMENT_H

#include <vector>
#include "StateSpace.hpp"
#include "TaxonSet.hpp"

// this class works like an interface
// it does not do any job
class SequenceAlignment {
  public:
    SequenceAlignment() {}

    SequenceAlignment(SequenceAlignment *from) {
        Ntaxa = from->Ntaxa;
        Nsite = from->Nsite;
        taxset = from->taxset;
        statespace = from->statespace;

        Data = from->Data;
    }

    double MissingFrac(int i) const {
        double n = 0;
        for (int k = 0; k < GetNtaxa(); k++) {
            if (Data[k][i] != unknown) {
                n++;
            }
        }
        return n / GetNtaxa();
    }

    SequenceAlignment(const int **inData, const std::string *names, int inNsite, const StateSpace *instatespace,
                      const TaxonSet *intaxset) {
        Nsite = inNsite;
        taxset = intaxset;
        Ntaxa = taxset->GetNtaxa();
        statespace = instatespace;

        Data.assign(Ntaxa,std::vector<int>(Nsite,0));
        for (int i = 0; i < Ntaxa; i++) {
            int mapi = taxset->GetTaxonIndex(names[i]);
            for (int j = 0; j < Nsite; j++) {
                Data[mapi][j] = inData[i][j];
            }
        }
    }

    virtual ~SequenceAlignment() = default;

    // the set of characters (A,C,G,T for nucleotides, etc..)
    const StateSpace *GetStateSpace() const { return statespace; }

    void Unclamp() {
        for (int i = 0; i < Ntaxa; i++) {
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = unknown;
            }
        }
    }

    int GetNstate() const { return statespace->GetNstate(); }

    // the list of taxa
    const TaxonSet *GetTaxonSet() const { return taxset; }

    int GetNsite() const { return Nsite; }

    int GetNtaxa() const { return taxset->GetNtaxa(); }

    bool isMissing(int taxon, int site) const { return Data[taxon][site] == -1; }

    bool AllMissingTaxon(int tax) const {
        bool ret = true;
        int site = 0;
        while ((site < GetNsite()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            site++;
        }
        return ret;
    }

    bool AllMissingTaxon(std::string taxname) const {
        int index = taxset->GetTaxonIndex(taxname);
        if (index == -1) {
            std::cerr << "error in all missing taxon: did not recognize " << taxname << '\n';
            exit(1);
        }
        return AllMissingTaxon(index);
    }

    bool AllMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            tax++;
        }
        return ret;
    }

    bool NoMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] != unknown);
            tax++;
        }
        return ret;
    }

    bool ConstantColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && (Data[tax][site] == unknown)) {
            tax++;
        }

        if (tax < GetNtaxa()) {
            int refstate = Data[tax][site];

            while ((tax < GetNtaxa()) && ret) {
                if (Data[tax][site] != -1) {
                    ret &= static_cast<int>(Data[tax][site] == refstate);
                }
                tax++;
            }
        }
        return ret;
    }

    void SetState(int taxon, int site, int state) { Data[taxon][site] = state; }
    int GetState(int taxon, int site) const { return Data[taxon][site]; }

    void ToStream(std::ostream &os) const;
    void ToStreamTriplet(std::ostream &os) const;
    int GetNonMissingTriplet() const;
    void ToFasta(std::ostream &os) const;

    // data fields

    int Ntaxa;
    int Nsite;
    const TaxonSet *taxset;
    const StateSpace *statespace;
    std::vector<std::vector<int> > Data;
};

class FileSequenceAlignment : public SequenceAlignment {
  public:
    FileSequenceAlignment(std::istream &is);
    FileSequenceAlignment(std::string filename);

  private:
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int ReadNexus(std::string filespec);
    int ReadSpecial(std::string filename);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);
};

#endif  // SEQUENCEALIGNMENT_H
