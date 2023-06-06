
#pragma once

#include "NodeArray.hpp"
#include "MPIBuffer.hpp"
#include "Random.hpp"

class PisArray : public SimpleNodeArray<double> {
  public:
    //! constructor 

    // specialized version
    // PisArray(int insize, const MultivariateBrownianTreeProcess& inX)

    // generic version (separation of concerns)
    PisArray(const Tree& intree, const NodeSelector<vector<double> >& inX, const double& invarfactor)
        : SimpleNodeArray<double>(intree), X(inX), varfactor(invarfactor), clamp(intree.GetNnode(),false)   {
        Sample();
    }

    ~PisArray() {}

    // indata: values of pi_s
    void ConditionOnData(const ContinuousData& indata)  {
    }

    const vector<double>& GetXLeafEntry(int i, int j) const {
        return X.GetVal(i).at(j);
    }

    double GetNl(intoo i) const   {
        double z = GetXLeafEntry(..,..);
        return ..;
    }

    double GetR(int i) const    {
    }

    double GetMean(int i) const {
        return ??;
    }

    double GetVar(int i) const {
        return ??;
    }

        // mean = shape / scale
        // var = shape / scale^2
        // shape = mean^2 / var
        // scale shape / mean

    double GetShape(int i) const { 
        double m = GetMean(i);
        double v = GetVar(i);
        return ??; 
    }

    double GetScale(int i) const { 
        return ??; 
    }

    //! sample all entries, given current shape and scale params
    void Sample() {
        for (int i = 0; i < GetNnode(); i++) {
            double shape = GetShape(i); ...
            (*this)[i] = Random::GammaSample(shape, scale);
        }
    }

    //! return total log prob (summed over the array) given current shape and
    //! scale params
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetNnode(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prob of one specific entry
    double GetLogProb(int index) const {
        // return non zero value only if pi_s is clamped (non-missing) for that species
        if (! clamp[index]) {
            return 0;
        }
        return Random::logGammaDensity(GetVal(index), shape, scale);
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetNnode(); i++) {
            m1 += GetVal(i);
        }
        m1 /= GetNnode();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetNnode(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetNnode();
        m2 /= GetNnode();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const NodeSelector<vector<double> >& X;
    const double& varfactor;
    vector<bool> clamp; // false / true
};
