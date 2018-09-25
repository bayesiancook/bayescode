#ifndef GENEMULTIGAMMA_H
#define GENEMULTIGAMMA_H

#include "Array.hpp"
#include "GammaSuffStat.hpp"

class MultiGamma : public SimpleArray<double>   {

    public:
    MultiGamma(const Selector<double>& inmean, const Selector<double>& ininvshape) :
        SimpleArray<double>(inmean.GetSize(),1.0), mean(inmean), invshape(ininvshape)   {
            Sample();
    }

    ~MultiGamma() {}

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            double shape = 1.0 / invshape.GetVal(i);
            double scale = shape / mean.GetVal(i);
            (*this)[i] = Random::GammaSample(shape, scale);
        }
    }

    //! return total log prob summed over all entries
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    //! return log prob for one entry
    double GetLogProb(int index) const {
        double shape = 1.0 / invshape.GetVal(index);
        double scale = shape / mean.GetVal(index);
        return Random::logGammaDensity(GetVal(index), shape, scale);
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
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const Selector<double> &mean;
    const Selector<double> &invshape;
};

class GeneIIDMultiGamma : public Array<MultiGamma>  {
    
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, shape parameter and center (frequency vector)
    GeneIIDMultiGamma(int inngene, const Selector<double>& inmean, const Selector<double>& ininvshape)
        : mean(inmean), invshape(ininvshape), size(inngene), array(inngene, (MultiGamma*) 0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new MultiGamma(mean, invshape);
        }
    }

    ~GeneIIDMultiGamma() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return size; }

    //! return total number of conditions
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access for the given gene
    const MultiGamma &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    MultiGamma &operator[](int gene) { return *array[gene]; }

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

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetSize(); gene++) {
            total += array[gene]->GetLogProb();
        }
        return total;
    }

    double GetLogProb(int gene, int cond) const {
        return array[gene]->GetLogProb(cond);
    }

  private:

    const Selector<double> &mean;
    const Selector<double> &invshape;
    int size;
    vector<MultiGamma*> array;
};

class MultiGammaSuffStat : public SimpleArray<GammaSuffStat>    {

    public:

    MultiGammaSuffStat(int insize) : SimpleArray<GammaSuffStat>(insize) {}
    ~MultiGammaSuffStat() {}

    void Clear()    {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i].Clear();
        }
    }

    void AddSuffStat(const GeneIIDMultiGamma& array)    {
        for (int gene=0; gene<array.GetSize(); gene++)  {
            for (int cond=0; cond<GetSize(); cond++)    {
                double x = array.GetVal(gene).GetVal(cond);
                (*this)[cond].AddSuffStat(x,log(x),1);
            }
        }
    }

    double GetLogProb(const Selector<double>& mean, const Selector<double>& invshape) const {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            double shape = 1.0 / invshape.GetVal(i);
            double scale = shape / mean.GetVal(i);
            total += GetVal(i).GetLogProb(shape,scale);
        }
        return total;
    }
};


#endif
