#pragma once

#include "Array.hpp"
#include "BranchArray.hpp"
#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"

class EpsRobustBranchIIDGamma : public SimpleBranchArray<double> {
  public:
    EpsRobustBranchIIDGamma(const Tree &intree, double inshape, double inscale, double inepsilon)
        : SimpleBranchArray<double>(intree), shape(inshape), scale(inscale), epsilon(inepsilon) {
        Sample();
    }

    ~EpsRobustBranchIIDGamma() {}

    double GetShape() const { return shape; }
    double GetScale() const { return scale; }
    double GetEpsilon() const { return epsilon; }

    void SetShape(double inshape) { shape = inshape; }
    void SetScale(double inscale) { scale = inscale; }
    void SetEpsilon(double inepsilon) { epsilon = inepsilon; }

    //! set all entries equal to inval
    void SetAllBranches(double inval) {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = inval;
        }
    }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            if (Random::Uniform() < epsilon)    {
                (*this)[i] = Random::CauchySample(1.0/scale);
            }
            else    {
                (*this)[i] = Random::GammaSample(shape, scale);
            }
        }
    }

    double Move(double tuning, int nrep, const BranchArray<PoissonSuffStat> &suffstatarray) {
        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<GetNbranch(); i++) {
                double logprob1 = GetLogProb(i) + suffstatarray.GetVal(i).GetLogProb(GetVal(i));
                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);
                double logh = m;
                (*this)[i] *= e;
                double logprob2 = GetLogProb(i) + suffstatarray.GetVal(i).GetLogProb(GetVal(i));
                double deltalogprob = logprob2 - logprob1 + logh;
                int accept = (log(Random::Uniform()) < deltalogprob);
                if (accept) {
                    nacc++;
                }
                else    {
                    (*this)[i] /= e;
                }
                ntot++;
            }
        }
        return nacc/ntot;
    }

    //! get total log prob summed over all branches
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! get log prob for a given branch
    double GetLogProb(int index) const { 
        double total = 0;
        total += log(1-epsilon) + Random::logGammaDensity(GetVal(index), shape, scale); 
        total += log(epsilon) + Random::logCauchyDensity(GetVal(index), 1.0/scale);
        return total;
    }

    //! get sum over all entries (name is rather specialized... could change..)
    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        return m1;
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetNbranch();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetNbranch();
        m2 /= GetNbranch();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    double shape;
    double scale;
    double epsilon;
};


