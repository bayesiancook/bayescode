
#ifndef PRODUCT_H
#define PRODUCT_H

/**
 * \brief The product of a Selector<double> and a double
 *
 * takes two arguments: an array l_j and a double multiplier a, and
 * returns, for index j, the product a * l_j.
 *
 */

class Product : public Selector<double> {
  public:
    //! \brief Constructor with branch array and multiplier
    Product(const Selector<double> &inbranchval, double inmulval)
        : branchval(inbranchval), mulval(inmulval), array(inbranchval.GetSize()) {}
    ~Product() {}

    virtual int GetSize() const override { return branchval.GetSize(); }

    virtual const double &GetVal(int index) const override {
        array[index] = mulval * branchval.GetVal(index);
        return array[index];
    }

    //! \brief specifies a new value for the multiplier
    void SetMulVal(double inmulval) { mulval = inmulval; }

  private:
    const Selector<double> &branchval;
    double mulval;
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

class ProductArray : public Array<Product> {
  public:
    //! Constructor: parameterized by the branch array and array
    ProductArray(const Selector<double> &inbranchval, const Selector<double> &insiteval)
        : branchval(inbranchval), siteval(insiteval), array(insiteval.GetSize(), (Product *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new Product(branchval, siteval.GetVal(gene));
        }
    }

    ~ProductArray() {
        for (int gene = 0; gene < GetSize(); gene++) { delete array[gene]; }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return siteval.GetSize(); }

    //! return total number of branches of the underlying tree
    int GetNcond() const { return array[0]->GetSize(); }

    //! const access to the Product for the given gene
    const Product &GetVal(int gene) const {
        // array[gene]->SetMulVal(siteval.GetVal(gene));
        return *array[gene];
    }

    //! non-const access to the Product for the given gene
    Product &operator[](int gene) { return *array[gene]; }

    //! global update of the arrays (needs to be called each time the array has
    //! changed)
    void Update() {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene]->SetMulVal(siteval.GetVal(gene));
        }
    }

  private:
    const Selector<double> &branchval;
    const Selector<double> &siteval;
    vector<Product *> array;
};

#endif
