#ifndef CODONSEQUENCEALIGNMENT_H
#define CODONSEQUENCEALIGNMENT_H

#include <cmath>
#include "CodonStateSpace.hpp"
#include "SequenceAlignment.hpp"

/**
 * \brief A codon sequence alignment (created from a nucleotide sequence
 * alignment)
 */

class CodonSequenceAlignment : public SequenceAlignment {
  public:
    //! \brief Constructor: takes a (const pointer to a) nucleotide sequence
    //! alignment and a genetic code.
    //!
    //! If force_stops is false, returns an error message when encountering stop
    //! codons -- otherwise, replace stop codons by missing entries. If any codon
    //! has any missing entry in the three nucleotide positions, then the whole
    //! codon is considered missing.
    CodonSequenceAlignment(
        SequenceAlignment *from, bool force_stops = false, GeneticCodeType type = Universal);

    ~CodonSequenceAlignment() /*override*/ = default;

    //! return the codon state space
    CodonStateSpace *GetCodonStateSpace() { return (CodonStateSpace *)(GetStateSpace()); }

    //! formatted output (simply, a nucleotide alignment)
    void ToStream(std::ostream &os);

  private:
    void ToStream(std::ostream &os, int pos);
    SequenceAlignment *DNAsource;
};

#endif
