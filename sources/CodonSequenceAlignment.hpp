#ifndef CODONSEQUENCEALIGNMENT_H
#define CODONSEQUENCEALIGNMENT_H

#include <cmath>
#include "CodonStateSpace.hpp"
#include "SequenceAlignment.hpp"

class CodonSequenceAlignment : public SequenceAlignment {
  public:
      /*
    CodonSequenceAlignment(CodonSequenceAlignment *from)
        : SequenceAlignment((SequenceAlignment *)from) {}
        */

    CodonSequenceAlignment(SequenceAlignment *from, bool force_stops = false,
                           GeneticCodeType type = Universal);

    ~CodonSequenceAlignment() /*override*/ = default;

    CodonStateSpace *GetCodonStateSpace() {
        // return static_cast<CodonStateSpace*>(statespace);
        return (CodonStateSpace *)(statespace);
    }

    void ToStream(std::ostream &os);
    void ToStream(std::ostream &os, int pos);

  private:
    SequenceAlignment *DNAsource;
};

#endif
