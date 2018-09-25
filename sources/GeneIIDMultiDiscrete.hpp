#ifndef GENEMULTIDISCRETE_H
#define GENEMULTIDISCRETE_H

#include "Array.hpp"

class MultiDiscrete : public SimpleArray<int>   {

    public:
    MultiDiscrete(const Selector<vector<double>>& inprob) : SimpleArray<int>(inprob.GetSize(),1), prob(inprob)  {
            Sample();
    }

    ~MultiDiscrete() {}

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::DrawFromDiscreteDistribution(prob.GetVal(i));
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
        return log(prob.GetVal(index)[GetVal(index)]);
    }

  protected:
    const Selector<vector<double>> &prob;
};

class GeneIIDMultiDiscrete : public Array<MultiDiscrete>  {
    
  public:
    //! constructor, parameterized by number of rows, of columns, dimension of the
    //! vectors, shape parameter and center (frequency vector)
    GeneIIDMultiDiscrete(int inngene, const Selector<vector<double>>& inprob) 
        : prob(inprob), size(inngene), array(inngene, (MultiDiscrete*) 0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new MultiDiscrete(prob);
        }
    }

    ~GeneIIDMultiDiscrete() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return size; }

    //! return total number of conditions
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access for the given gene
    const MultiDiscrete &GetVal(int gene) const { return *array[gene]; }

    //! non-const access for the given gene
    MultiDiscrete &operator[](int gene) { return *array[gene]; }

    //! return total log prob (over all genes and over all branches)
    double GetLogProb() const {
        double total = 0;
        for (int gene = 0; gene < GetSize(); gene++) {
            total += array[gene]->GetLogProb();
        }
        return total;
    }

  private:

    const Selector<vector<double>> &prob;
    int size;
    vector<MultiDiscrete*> array;
};

#endif
