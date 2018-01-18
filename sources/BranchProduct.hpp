
#ifndef BRANCHPRODUCT_H
#define BRANCHPRODUCT_H

class BranchProduct : public BranchSelector<double>  {

    public:

    BranchProduct(const BranchSelector<double>& inbranchval, double inmulval) : branchval(inbranchval), mulval(inmulval), array(inbranchval.GetNbranch()) {}
    ~BranchProduct() {}

	virtual const Tree& GetTree() const {return branchval.GetTree();}
	virtual const double& GetVal(int index) const {
        array[index] = mulval * branchval.GetVal(index);
        return array[index];
    }

    void SetMulVal(double inmulval)   {
        mulval = inmulval;
    }

    private:

    const BranchSelector<double>& branchval;
    double mulval;
    mutable vector<double> array;
};

class BranchProductArray : public Array<BranchProduct>    {

    public:

    //! constructor: parameterized by the number of genes, the tree, the means over branches and the shape parameter
    BranchProductArray(const BranchSelector<double>& inbranchval, const Selector<double>& insiteval) : branchval(inbranchval), siteval(insiteval), array(insiteval.GetSize(), (BranchProduct*) 0) {
        for (int gene=0; gene<GetSize(); gene++)    {
            array[gene] = new BranchProduct(branchval,siteval.GetVal(gene));
        }
    }

    ~BranchProductArray()  {
        for (int gene=0; gene<GetSize(); gene++)    {
            delete array[gene];
        }
    }

    //! return total number of entries (number of genes)
    int GetSize() const {
        return siteval.GetSize();
    }

    //! return total number of branches of the underlying tree
    int GetNbranch() const {
        return array[0]->GetNbranch();
    }

    //! const access to the BranchProduct for the given gene
    const BranchProduct& GetVal(int gene) const {
        // array[gene]->SetMulVal(siteval.GetVal(gene));
        return *array[gene];
    }

    //! non-const access to the BranchProduct for the given gene
    BranchProduct& operator[](int gene)  {
        return *array[gene];
    }

    void Update()   {
        for (int gene=0; gene<GetSize(); gene++)    {
            array[gene]->SetMulVal(siteval.GetVal(gene));
        }
    }

    private:

    const BranchSelector<double>& branchval;
    const Selector<double>& siteval;
    vector<BranchProduct*> array;
};

/*
class BranchSiteProduct : public BranchSiteSelector<double> {

    public:

    BranchSiteProduct(const BranchSelector<double>& inbranchval, const Selector<double>& insiteval) : branchval(inbranchval), siteval(insiteval), array(insiteval.GetSize(),vector<double>(inbranchval.GetNbranch(),1.0)) {}

    ~BranchSiteProduct() {}

	virtual const Tree& GetTree() const  {return branchval.GetTree();}
	virtual int GetSize() const {return siteval.GetSize();}

	virtual const double& GetVal(int branch, int site) const {
        array[site][branch] = siteval.GetVal(site) * branchval.GetVal(branch);
        return array[site][branch];
    }

    private:

    const BranchSelector<double>& branchval;
    const Selector<double>& siteval;
    mutable vector<vector<double> > array;
};
*/

#endif

