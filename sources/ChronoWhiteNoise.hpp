#pragma once

#include "BranchArray.hpp"
#include "PoissonSuffStat.hpp"

/**
 * \brief A tree-structured branch-wise array of Gamma variables, with
 * branch-specific means but same shape parameter
 *
 * One should be careful about the fact that the shape parameter is given by
 * copy (not by ref) to the array. Thus, each time the shape parameter is
 * modified during the MCMC, the new value should be given to the array (using
 * the SetShape method).
 */
class ChronoGammaWhiteNoise : public SimpleBranchArray<double> {
  public:

    // mode 0   : ugam
    // mode 1/2 : wn of mean 1 and variance per unit of time = 1/shape
    // mode 1   : realized branch variable: integral of wn over time interval = blmean
    // mode 2   : realized branch variable: mean     of wn over time interval = blmean

    ChronoGammaWhiteNoise(const Tree &intree, const NodeSelector<double> &inchrono, double inshape, int inmode = 2)
        : SimpleBranchArray<double>(intree), chrono(inchrono), shape(inshape) {
        mode = inmode;
        Sample();
    }

    ~ChronoGammaWhiteNoise() {}

    void SetMode(int inmode)    {
        mode = inmode;
    }

    double GetShape() const { return shape; }

    void SetShape(double inshape) { shape = inshape; }

    double GetAlpha(const Link* from) const {
        if (mode == 3)  {
            return shape;
        }
        if (mode == 2)  {
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            return shape * dt;
        }
        if (mode == 1)   {
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            return shape * dt;
        }
        return shape;
    }

    double GetBeta(const Link* from) const {
        if (mode == 3)  {
            return shape;
        }
        if (mode == 2)  {
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            return shape * dt;
        }
        if (mode == 1)   {
            return shape;
        }
        double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
        return shape / dt;
    }
        
    //! sample all entries from prior
    void Sample() {
        RecursiveSample(GetTree().GetRoot());
    }

    void RecursiveSample(const Link* from)  {
        if (! from->isRoot())   {
            (*this)[from->GetBranch()->GetIndex()] = Random::GammaSample(GetAlpha(from), GetBeta(from));
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveSample(link->Out());
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void GibbsResample(const PoissonSuffStatBranchArray &suffstatarray) {
        RecursiveGibbsResample(GetTree().GetRoot(), suffstatarray);
    }

    void RecursiveGibbsResample(const Link* from, const PoissonSuffStatBranchArray &suffstatarray) {
        if (! from->isRoot())   {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(from->GetBranch()->GetIndex());
            double tmp = Random::GammaSample(GetAlpha(from) + suffstat.GetCount(), GetBeta(from) + suffstat.GetBeta());
            if (! tmp)  {
                cerr << "null sample in white noise: " << GetAlpha(from) << '\t' << GetBeta(from) << '\t' << suffstat.GetCount() << '\t' << suffstat.GetBeta() << '\n';
                cerr << mode << '\t' << shape << '\n';
                exit(1);
            }
            (*this)[from->GetBranch()->GetIndex()] = tmp;
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveGibbsResample(link->Out(), suffstatarray);
        }
    }

    //! return total log prob summed over all entries
    double GetLogProb() {
        return RecursiveGetLogProb(GetTree().GetRoot());
    }

    double RecursiveGetLogProb(const Link* from)    {
        double total = 0;
        if (! from->isRoot())   {
            total += GetBranchLogProb(from);
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += RecursiveGetLogProb(link->Out());
        }
        return total;
    }

    //! return log prob for one entry
    double GetBranchLogProb(const Link* from) const {
        return Random::logGammaDensity(GetVal(from->GetBranch()->GetIndex()), GetAlpha(from), GetBeta(from));
    }

    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1 / GetNbranch();
    }

  protected:
    const NodeSelector<double> &chrono;
    double shape;
    int mode;
};

