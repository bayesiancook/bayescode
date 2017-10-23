#ifndef SBP_H
#define SBP_H

#include "Array.hpp"

class StickBreakingProcess : public SimpleArray<double> {

    public:

    StickBreakingProcess(int inncat, double inkappa) : SimpleArray<double>(inncat), kappa(inkappa), V(inncat)    {
        Sample();
    }

    ~StickBreakingProcess() {}

    void SetKappa(double inkappa)   {
        kappa = inkappa;
    }

    const vector<double>& GetBetaVariates() const {
        return V;
    }

    void SwapComponents(int cat1, int cat2)   {
        Swap(cat1,cat2);
        double tmp = V[cat1];
        V[cat1] = V[cat2];
        V[cat2] = tmp;
    }

    void Sample()   {
        double cumulProduct = 1.0;
        double totweight = 0;
        for (int k=0; k<GetSize(); k++)	{
            double x = Random::sGamma(1.0);
            double y = Random::sGamma(kappa);
            double v = x / (x+y);
            V[k] = v;
            if (k == GetSize() - 1)	{
                V[k] = 1;
                v = 1;
            }
            (*this)[k] = v * cumulProduct;
            cumulProduct *= (1 - v);	
            totweight += (*this)[k];
        }
    }


    void GibbsResample(const vector<int>& occupancy)	{

        int remainingOcc = 0;
        for (unsigned int i=0; i<occupancy.size(); i++) {
            remainingOcc += occupancy[i];
        }

        double cumulProduct = 1.0;
        double totweight = 0;
        for (int k=0; k<GetSize(); k++)	{
            remainingOcc -= occupancy[k];
            double x = Random::sGamma(1 + occupancy[k]);
            double y = Random::sGamma(kappa + remainingOcc);
            double v = x / (x+y);
            V[k] = v;
            if (k == GetSize() - 1)	{
                V[k] = 1;
                v = 1;
            }
            (*this)[k] = v * cumulProduct;
            cumulProduct *= (1 - v);
            totweight += (*this)[k];
        }
    }

    double GetLogProb() const   {
        double total = 0;
        for (int k=0; k<GetSize()-1; k++)	{
            total += Random::logBetaDensity(V[k],1.0,kappa);
        }
        return total;
    }

    double GetLogProb(double kappa) const   {
        double total = 0;
        for (int k=0; k<GetSize()-1; k++)	{
            total += Random::logBetaDensity(V[k],1.0,kappa);
        }
        return total;
    }

    double GetMarginalLogProb(const vector<int>& occupancy) const   {

        int remainingOcc = 0;
        for (unsigned int i=0; i<occupancy.size(); i++) {
            remainingOcc += occupancy[i];
        }

        double total = 0;
        for (int k=0; k<GetSize(); k++)	{
            if (remainingOcc)	{
                remainingOcc -= occupancy[k];
                total += log(kappa) + Random::logGamma(1 + occupancy[k]) + Random::logGamma(kappa + remainingOcc) - Random::logGamma(1 + kappa + occupancy[k] + remainingOcc);
            }
        }
        if (remainingOcc)	{
            cerr << "error in allocation count\n";
            exit(1);
        }
        return total;
    }

    private:

    double kappa;
    vector<double> V;

};

#endif

