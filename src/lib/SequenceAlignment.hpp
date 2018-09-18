#ifndef SEQUENCEALIGNMENT_H
#define SEQUENCEALIGNMENT_H

#include <vector>
#include "StateSpace.hpp"
#include "TaxonSet.hpp"

/**
 * \brief Generic interface for a multiple sequence alignment
 */

class SequenceAlignment {
  public:
    //! default constructor
    SequenceAlignment()
        : Ntaxa(0),
          Nsite(0),
          taxset(nullptr),
          statespace(nullptr),
          owntaxset(true),
          ownstatespace(true) {}

    //! copy constructor
    SequenceAlignment(const SequenceAlignment &from)
        : Ntaxa(from.Ntaxa),
          Nsite(from.Nsite),
          taxset(from.taxset),
          statespace(from.statespace),
          owntaxset(false),
          ownstatespace(false),
          Data(from.Data) {}

    virtual ~SequenceAlignment() {
        if (owntaxset) {
            delete taxset;
            taxset = nullptr;
        }
        if (ownstatespace) {
            delete statespace;
            statespace = nullptr;
        }
    }

    //! return the state space (A,C,G,T for nucleotides, etc..)
    const StateSpace *GetStateSpace() const { return statespace; }

    //! return size of state space
    int GetNstate() const { return statespace->GetNstate(); }

    //! return the set of taxa
    const TaxonSet *GetTaxonSet() const { return taxset; }

    //! return the number of aligned positions
    int GetNsite() const { return Nsite; }

    //! return the number of aligned positions such as printed out (Nsite for
    //! nucleotide models, 3*Nsite for codon models)
    int GetPrintNsite() const { return GetStateSpace()->GetSymbolLength() * Nsite; }

    //! return the number of taxa (number of aligned sequences)
    int GetNtaxa() const { return taxset->GetNtaxa(); }

    // return state for this taxon at that site (return -1 if missing entry)
    int GetState(int taxon, int site) const { return Data[taxon][site]; }

    //! whether or not entry is missing for this taxon at that site
    bool isMissing(int taxon, int site) const { return Data[taxon][site] == -1; }

    //! Phylip-like formatted output to stream
    void ToStream(std::ostream &os) const;

    //! set the state to a new value (note: should really re-consider this option,
    //! currently used by PhyloProcess to simulate new data)
    void SetState(int taxon, int site, int state) { Data[taxon][site] = state; }

    //! return empirical frequencies into a vector
    std::vector<double> GetEmpiricalFreq() const;

  protected:
    bool AllMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            tax++;
        }
        return ret;
    }

  private:
    // replace all entries by missing entries
    void Unclamp() {
        for (int i = 0; i < Ntaxa; i++) {
            for (int j = 0; j < Nsite; j++) { Data[i][j] = unknown; }
        }
    }

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
        while ((tax < GetNtaxa()) && (Data[tax][site] == unknown)) { tax++; }

        if (tax < GetNtaxa()) {
            int refstate = Data[tax][site];

            while ((tax < GetNtaxa()) && ret) {
                if (Data[tax][site] != -1) { ret &= static_cast<int>(Data[tax][site] == refstate); }
                tax++;
            }
        }
        return ret;
    }

    void ToStreamTriplet(std::ostream &os) const;
    int GetNonMissingTriplet() const;
    void ToFasta(std::ostream &os) const;

    // data fields

  protected:
    int Ntaxa;
    int Nsite;
    const TaxonSet *taxset;
    const StateSpace *statespace;
    bool owntaxset;
    bool ownstatespace;
    std::vector<std::vector<int>> Data;
};

/**
 * \brief A sequence alignment created by reading from a file (Phylip-like or
 * Nexus format)
 */

class FileSequenceAlignment : public SequenceAlignment {
  public:
    FileSequenceAlignment(std::string filename);
    FileSequenceAlignment(std::istream &is);

  private:
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int ReadNexus(std::string filespec);
    int ReadSpecial(std::string filename);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    void ReadPhylipSequentialFromStream(std::istream &is);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);
};

#endif  // SEQUENCEALIGNMENT_H
