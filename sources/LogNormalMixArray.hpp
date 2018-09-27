
#ifndef LOGNORMALMIXARRAY_H
#define LOGNORMALMIXARRAY_H

#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Product.hpp"
#include "Random.hpp"
#include "GeneIIDMultiGamma.hpp"
#include "GeneIIDMultiDiscrete.hpp"


class LogNormalMixArray : public SimpleArray<double> {
  public:
    LogNormalMixArray(const Selector<double> &inmean, const Selector<int>& inalloc, const Selector<double>& indevpos, const Selector<double>& indevneg) : SimpleArray<double>(inmean.GetSize()), mean(inmean), alloc(inalloc), devpos(indevpos), devneg(indevneg)    {
        Update();
    }

    ~LogNormalMixArray() {}

    //! sample all entries from prior
    void Update() {
        for (int i = 0; i < GetSize(); i++) {
            double x = log(mean.GetVal(i));
            if (alloc.GetVal(i) == 0)   {
                x -= devneg.GetVal(i);
            }
            else if (alloc.GetVal(i) == 2)  {
                x += devpos.GetVal(i);
            }
            (*this)[i] = exp(x);
        }
    }

    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
        }
        return m1 / GetSize();
    }

    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            double tmp = GetVal(i);
            m1 += tmp;
            m2 += tmp * tmp;
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

    double GetMeanLog() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += log(GetVal(i));
        }
        return m1 / GetSize();
    }

    double GetVarLog() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            double tmp = log(GetVal(i));
            m1 += tmp;
            m2 += tmp * tmp;
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const Selector<double>& mean;
    const Selector<int>& alloc;
    const Selector<double>& devpos;
    const Selector<double>& devneg;
};

class GeneIIDLogNormalMixArray : public Array<LogNormalMixArray> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    GeneIIDLogNormalMixArray(const ProductArray &inmean, const GeneIIDMultiDiscrete& inalloc, const GeneIIDMultiGamma& indevpos, const GeneIIDMultiGamma& indevneg) : mean(inmean), alloc(inalloc), devpos(indevpos), devneg(indevneg), array(inmean.GetSize(), (LogNormalMixArray *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new LogNormalMixArray(mean.GetVal(gene), alloc.GetVal(gene), devpos.GetVal(gene), devneg.GetVal(gene));
        }
    }

    ~GeneIIDLogNormalMixArray() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
        }
    }

    void Update()   {
        for (int gene=0; gene<GetSize(); gene++)   {
            array[gene]->Update();
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return mean.GetSize(); }

    //! return total number of conditions
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access for the given gene
    const LogNormalMixArray &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    LogNormalMixArray &operator[](int gene) { return *array[gene]; }

    //! return mean tree length over genes
    double GetMean() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = array[g]->GetVal(j);
                mean += tmp;
            }
            mean /= GetSize();
            tot += mean;
        }
        return tot;
    }

    //! return variance of tree length across genes
    double GetVar() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = array[g]->GetVal(j);
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= GetSize();
            var /= GetSize();
            var -= mean * mean;
            tot += var;
        }
        tot /= GetSize();
        return tot;
    }

    //! return mean tree length over genes
    double GetMeanLog() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = log(array[g]->GetVal(j));
                mean += tmp;
            }
            mean /= GetSize();
            tot += mean;
        }
        return tot;
    }

    //! return variance of tree length across genes
    double GetVarLog() const {
        double tot = 0;
        for (int j = 0; j < GetNcond(); j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < GetSize(); g++) {
                double tmp = log(array[g]->GetVal(j));
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= GetSize();
            var /= GetSize();
            var -= mean * mean;
            tot += var;
        }
        tot /= GetSize();
        return tot;
    }

  private:
    const ProductArray &mean;
    const GeneIIDMultiDiscrete& alloc;
    const GeneIIDMultiGamma& devpos;
    const GeneIIDMultiGamma& devneg;
    vector<LogNormalMixArray *> array;
};

#endif
