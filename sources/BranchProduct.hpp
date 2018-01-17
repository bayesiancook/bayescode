
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

/*
class BranchSiteProduct : public BranchSiteSelector<double> {

    public:

    BranchSiteProduct(const BranchSelector<double>& inbranchval, const Selector<double>& insiteval) : branchval(inbranchval), siteval(insiteval) {}
    ~BranchSiteProduct() {}

	virtual const Tree& GetTree() const  {return branchval.GetTree();}
	virtual int GetSize() const {return siteval.GetSize();}
	virtual const double GetVal(int branch, int site) const {return siteval.GetVal(site) * branchval.GetVal(branch);}

    private:

    const BranchSelector<double>& branchval;
    const Selector<double>& siteval;

};
*/

#endif

