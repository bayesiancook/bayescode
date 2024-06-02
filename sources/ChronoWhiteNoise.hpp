#pragma once

#include "BranchArray.hpp"
#include "PoissonSuffStat.hpp"

/**
 * \brief A tree-structured branch-wise array of Gamma variables, with
 * branch-specific means but same variance parameter
 *
 * One should be careful about the fact that the var parameter is given by
 * copy (not by ref) to the array. Thus, each time the var parameter is
 * modified during the MCMC, the new value should be given to the array (using
 * the SetShape method).
 */
class ChronoGammaWhiteNoise : public SimpleBranchArray<double> {
  public:

    // mode 1   : variance var1
    // mode 2   : variance var1*dt
    // mode 3   : variance var1/dt
    // mode 4   : variance var1/dt + var2 + var3*dt

    ChronoGammaWhiteNoise(const Tree &intree, const NodeSelector<double> &inchrono, int inmode)
        : SimpleBranchArray<double>(intree), chrono(inchrono), var1(1.0), var2(1.0), var3(1.0)  {
        mode = inmode;
        Sample();
    }

    ~ChronoGammaWhiteNoise() {}

    void SetVar(double invar) { var1 = invar; }

    void SetVar(double in1, double in2, double in3) {
        var1 = in1;
        var2 = in2;
        var3 = in3;
    }

    double GetAlpha(const Link* from) const {
        double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
        switch(mode)    {
            case 1:
                return 1.0/var1;
                break;
            case 2:
                return 1.0 / (var1/dt);
                break;
            case 3:
                return 1.0 / (var1*dt);
                break;
            case 4:
                return 1.0 / (var1/dt + var2 + var3*dt);
                break;
            default:
                cerr << "unrecognized mode in white noise\n";
                exit(1);
        }
    }

    //! sample all entries from prior
    void Sample() {
        RecursiveSample(GetTree().GetRoot());
    }

    void RecursiveSample(const Link* from)  {
        if (! from->isRoot())   {
            (*this)[from->GetBranch()->GetIndex()] = Random::GammaSample(GetAlpha(from), GetAlpha(from));
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveSample(link->Out());
        }
    }

    //! resample entries based on a BranchArray of PoissonSuffStat
    void GibbsResample(const PoissonSuffStatBranchArray &suffstatarray, double b=1) {
        RecursiveGibbsResample(GetTree().GetRoot(), suffstatarray, b);
    }

    void RecursiveGibbsResample(const Link* from, const PoissonSuffStatBranchArray &suffstatarray, double b=1) {
        if (! from->isRoot())   {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(from->GetBranch()->GetIndex());
            double tmp = Random::GammaSample(GetAlpha(from) + b*suffstat.GetCount(), GetAlpha(from) + b*suffstat.GetBeta());
            if (! tmp)  {
                tmp = (GetAlpha(from) + b*suffstat.GetCount()) / (GetAlpha(from) + b*suffstat.GetBeta());
                /*
                cerr << "null sample in white noise: " << GetAlpha(from) << '\t' << suffstat.GetCount() << '\t' << suffstat.GetBeta() << '\n';
                cerr << (GetAlpha(from) + suffstat.GetCount()) / (GetAlpha(from) + suffstat.GetBeta()) << '\n';
                exit(1);
                */
            }
            (*this)[from->GetBranch()->GetIndex()] = tmp;
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveGibbsResample(link->Out(), suffstatarray);
        }
    }

    void VariancePartition(double& v1, double& v2, double& v3)  {
        RecursiveVariancePartition(GetTree().GetRoot(), v1, v2, v3);
    }

    void RecursiveVariancePartition(const Link* from, double& v1, double& v2, double& v3)   {
        if (! from->isRoot())   {
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            double tmp1 = var1/dt;
            double tmp2 = var2;
            double tmp3 = var3*dt;
            double tot = tmp1 + tmp2 + tmp3;
            v1 += tmp1/tot;
            v2 += tmp2/tot;
            v3 += tmp3/tot;
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveVariancePartition(link->Out(), v1, v2, v3);
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
        return Random::logGammaDensity(GetVal(from->GetBranch()->GetIndex()), GetAlpha(from), GetAlpha(from));
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
    double var1, var2, var3;
    int mode;
};

