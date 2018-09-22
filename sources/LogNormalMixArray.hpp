
#ifndef LOGNORMALMIXARRAY_H
#define LOGNORMALMIXARRAY_H

#include "MPIBuffer.hpp"
#include "PoissonSuffStat.hpp"
#include "Sum.hpp"
#include "Random.hpp"

class LogNormalMixArray : public SimpleArray<double> {
  public:
    LogNormalMixArray(const Selector<double> &inmean, const Selector<double>& inpipos, const Selector<double>& inmeanpos, const Selector<double>& ininvshapepos, const Selector<double>& inpineg, const Selector<double>& inmeanneg, const Selector<double>& ininvshapeneg) 
        : SimpleArray<double>(inmean.GetSize()), mean(mean), pipos(inpipos), meanpos(inmeanpos), invshapepos(ininvshapepos), pineg(inpineg), meanneg(inmeanneg), invshapeneg(ininvshapeneg) {
        Sample();
    }

    ~LogNormalMixArray() {}

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            double logv = 0;
            double pi = pipos.GetVal(i) + pineg.GetVal(i);
            if (pi > 1.0)   {
                cerr << "error in LogNormalMixArray::Sample: prob > 1\n";
                exit(1);
            }
            double u = Random::Uniform();
            logv = mean.GetVal(i);

            if (u < pipos.GetVal(i))   {
                double shape = 1.0 / invshapepos.GetVal(i);
                double scale = shape / meanpos.GetVal(i);
                logv += Random::GammaSample(shape, scale);
            }
            else if (u < pi)    {
                double shape = 1.0 / invshapeneg.GetVal(i);
                double scale = shape / meanneg.GetVal(i);
                logv -= Random::GammaSample(shape, scale);
            }
            (*this)[i] = exp(logv);
        }
    }

    //! resample entries based on a Array of PoissonSuffStat
    void MultipleTryResample(const Selector<PoissonSuffStat> &suffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
            cerr << "do something here\n";
            exit(1);
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
        double dev = log(GetVal(index)) - mean.GetVal(index);
        double ret = 0;
        if (dev > 0)    {
            ret += log(pipos.GetVal(index));
            double shape = 1.0 / invshapepos.GetVal(index);
            double scale = shape / meanpos.GetVal(index);
            ret += Random::logGammaDensity(dev, shape, scale);
        }
        else if (dev < 0)   {
            ret += log(pineg.GetVal(index));
            double shape = 1.0 / invshapeneg.GetVal(index);
            double scale = shape / meanneg.GetVal(index);
            ret += Random::logGammaDensity(-dev, shape, scale);
        }
        else    {
            double pi = pipos.GetVal(index) + pineg.GetVal(index);
            if (pi >= 1.0)  {
                cerr << "error in LogNormalMixArray:GetLogProb: pi > 1.0\n";
                exit(1);
            }
            ret += log(1.0 - pi);
        }
        return ret;
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
    const Selector<double>& pipos;
    const Selector<double>& meanpos;
    const Selector<double>& invshapepos;
    const Selector<double>& pineg;
    const Selector<double>& meanneg;
    const Selector<double>& invshapeneg;
};

class LogNormalMixBidimArray : public Array<LogNormalMixArray> {
  public:
    //! constructor: parameterized by the number of genes, the tree, the means
    //! over branches and the shape parameter
    LogNormalMixBidimArray(const SumArray &inmean, const Selector<double>& inpipos, const Selector<double>& inmeanpos, const Selector<double>& ininvshapepos, const Selector<double>& inpineg, const Selector<double>& inmeanneg, const Selector<double>& ininvshapeneg)
        : mean(inmean), pipos(inpipos), meanpos(inmeanpos), invshapepos(ininvshapepos), pineg(inpineg), meanneg(inmeanneg), invshapeneg(ininvshapeneg),
          array(inmean.GetSize(), (LogNormalMixArray *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new LogNormalMixArray(mean.GetVal(gene), pipos, meanpos, invshapepos, pineg, meanneg, invshapeneg);
        }
    }

    ~LogNormalMixBidimArray() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
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

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetSize(); gene++) {
            total += array[gene]->GetLogProb();
        }
        return total;
    }

  private:
    const SumArray &mean;
    const Selector<double>& pipos;
    const Selector<double>& meanpos;
    const Selector<double>& invshapepos;
    const Selector<double>& pineg;
    const Selector<double>& meanneg;
    const Selector<double>& invshapeneg;
    vector<LogNormalMixArray *> array;
};

#endif
