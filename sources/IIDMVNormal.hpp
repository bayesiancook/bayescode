
#ifndef IIDMVNORMAL_H
#define IIDMVNORMAL_H

// #include "Array.hpp"
#include "BidimArray.hpp"

class BidimIIDMVNormal : public SimpleBidimArray<vector<double> >    {

    public:

    BidimIIDMVNormal(int inncol, int indim, const ConstArray<double>& invar) : 
        SimpleBidimArray(invar.GetSize(),inncol,vector<double>(indim,0)), dim(indim), var(invar) {
        Sample();
    }

    int GetDim() const {return dim;}

    void Sample()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                Sample(i,j);
            }
        }
    }

    void Sample(int i, int j)   {
        vector<double>& x = (*this)(i,j);
        double v = var.GetVal(i);
        for (int k=0; k<GetDim(); k++) {
            x[k] = sqrt(v)*Random::sNormal();
        }
    }

    double GetMeanVar(int k) const  {

        double mean = 0;
        for (int j=0; j<GetNcol(); j++) {
            mean += GetVar(k,j);
        }
        mean /= GetNcol();
        return mean;
    }

    double GetVar(int k, int j) const   {

        double mean = 0;
        double var = 0;
        const vector<double>& x = GetVal(k,j);
        for (int l=0; l<GetDim(); l++)  {
            mean += x[l];
            var += x[l]*x[l];
        }
        mean /= GetDim();
        var /= GetDim();
        var -= mean*mean;
        return var;
    }

    double GetLogProb() const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetColumnLogProb(j);
        }
        return total;
    }

    double GetRowLogProb(int i) const {
        double total = 0;
        for (int j=0; j<GetNcol(); j++)   {
            total += GetLogProb(i,j);
        }
        return total;
    }

    double GetColumnLogProb(int j) const  {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            total += GetLogProb(i,j);
        }
        return total;
    }

    double GetColumnLogProb(int j, const vector<int>& flag) const   {
        double total = 0;
        for (int i=0; i<GetNrow(); i++)  {
            if (flag[i])    {
                total += GetLogProb(i,j);
            }
        }
        return total;
    }

    double GetLogProb(int i, int j) const {
        double total = 0;
        const vector<double>& x = GetVal(i,j);
        double v = var.GetVal(i);
        for (int k=0; k<GetDim(); k++) {
            total += x[k]*x[k];
        }
        return -0.5*(total/v + GetDim()*log(v));
    }

    protected:

    int dim;
    const ConstArray<double>& var;
};


/*
class IIDMVNormal : public SimpleArray<vector<double >    {

    public:

    IIDMVNormal(int insize, int indim, double invar) : SimpleArray<vector<double> >(insize,vector<double>(indim,0)), dim(indim), var(invar)  {
        Sample();
    }
    
    ~IIDMVNormal() {}

    void Sample()   {
        for (int i=0; i<GetSize(); i++) {
            for (int j=0; j<GetDim(); j++) {
                (*this)[i][j] = sqrt(var)*Random::sNormal();
            }
        }
    }

    double GetLogProb() {
        double total = 0;
        for (int i=0; i<GetSize(); i++) {
            total += GetLogProb(i);
        }
        return total;
    }

    double GetLogProb(int i)    {
        double total = 0;
        for (int j=0; j<GetDim(); j++)  {
            double tmp = GetVal(i)[j];
            total += tmp*tmp;
        }
        double ret = -0.5*(total/var + GetDim()*log(var));
    }

    int GetDim()    {
        return dim;
    }

    void SetVar(double invar)   {
        var = invar;
    }

    protected:

    int dim;

    double var;

};
*/

#endif

