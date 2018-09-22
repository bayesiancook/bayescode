
#ifndef SUM_H
#define SUM_H

/**
 * \brief The product of a Selector<double> and a double
 *
 * takes two arguments: an array l_j and a double multiplier a, and
 * returns, for index j, the product a * l_j.
 *
 */

class Sum : public Selector<double> {
  public:
    //! \brief Constructor with branch array and multiplier
    Sum(const Selector<double> &inbranchval, double ingeneval)
        : branchval(inbranchval), geneval(ingeneval), array(inbranchval.GetSize()) {}
    ~Sum() {}

    virtual int GetSize() const override { return branchval.GetSize(); }

    virtual const double &GetVal(int index) const override {
        array[index] = geneval + branchval.GetVal(index);
        return array[index];
    }

    //! \brief specifies a new value for the multiplier
    void SetGeneVal(double ingeneval) { geneval = ingeneval; }

  private:
    const Selector<double> &branchval;
    double geneval;
    mutable vector<double> array;
};

/**
 * \brief The product of a Selector<double> and a Selector<double>
 *
 * Takes two arguments: a Array l_j and an Array r_i, and
 * returns, for branch j, and item i, the product r_i * l_j.
 *
 * Used in OmegaModel: the mean of omega_gj, for gene g and branch j, is equal
 * to w_g * branchv_j
 */

class SumArray : public Array<Sum> {
  public:
    //! Constructor: parameterized by the branch array and array
    SumArray(const Selector<double> &inbranchval, const Selector<double> &ingeneval)
        : branchval(inbranchval), geneval(ingeneval), array(ingeneval.GetSize(), (Sum *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new Sum(branchval, geneval.GetVal(gene));
        }
    }

    ~SumArray() {
        for (int gene = 0; gene < GetSize(); gene++) {
            delete array[gene];
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return geneval.GetSize(); }

    //! return total number of branches of the underlying tree
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access to the Sum for the given gene
    const Sum &GetVal(int gene) const {
        // array[gene]->SetGeneVal(geneval.GetVal(gene));
        return *array[gene];
    }

    //! non-const access to the Sum for the given gene
    Sum &operator[](int gene) { return *array[gene]; }

    //! global update of the arrays (needs to be called each time the array has
    //! changed)
    void Update() {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene]->SetGeneVal(geneval.GetVal(gene));
        }
    }

  private:
    const Selector<double> &branchval;
    const Selector<double> &geneval;
    vector<Sum *> array;
};

#endif
