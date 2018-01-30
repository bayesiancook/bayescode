#ifndef SBP_H
#define SBP_H

#include "Array.hpp"
#include "OccupancySuffStat.hpp"
#include "Permutation.hpp"

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

    void GibbsResample(const OccupancySuffStat& occupancy)  {

        int remainingOcc = 0;
        for (int i=0; i<occupancy.GetSize(); i++) {
            remainingOcc += occupancy.GetVal(i);
        }

        double cumulProduct = 1.0;
        double totweight = 0;
        for (int k=0; k<GetSize(); k++)	{
            remainingOcc -= occupancy.GetVal(k);
            double x = Random::sGamma(1 + occupancy.GetVal(k));
            double y = Random::sGamma(kappa + remainingOcc);
            double v = x / (x+y);
            if (! v)    {
                v = 1e-50;
            }
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

    double GetMarginalLogProb(const OccupancySuffStat& occupancy) const {

        int remainingOcc = 0;
        for (int i=0; i<occupancy.GetSize(); i++) {
            remainingOcc += occupancy.GetVal(i);
        }

        double total = 0;
        for (int k=0; k<GetSize(); k++)	{
            if (remainingOcc)	{
                remainingOcc -= occupancy.GetVal(k);
                total += log(kappa) + Random::logGamma(1 + occupancy.GetVal(k)) + Random::logGamma(kappa + remainingOcc) - Random::logGamma(1 + kappa + occupancy.GetVal(k) + remainingOcc);
            }
        }
        if (remainingOcc)	{
            cerr << "error in allocation count\n";
            exit(1);
        }
        return total;
    }

    void LabelSwitchingMove(int nrep, OccupancySuffStat& occupancy, Permutation& permut)   {
        MoveOccupiedCompAlloc(nrep,occupancy,permut);
        MoveAdjacentCompAlloc(nrep,occupancy,permut);
    }

    double MoveOccupiedCompAlloc(int k0, OccupancySuffStat& occupancy, Permutation& permut)	{

        int nrep = (int) (k0 * kappa);
        GibbsResample(occupancy);
        double total = 0.0;
        int Nocc = occupancy.GetNcluster();
        if (Nocc != 1)	{
            for (int i=0; i<nrep; i++)	{
                int occupiedComponentIndices[Nocc];
                int j=0;
                for (int k=0; k<GetSize(); k++)	{
                    if (occupancy[k] != 0)	{
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc)	{
                    cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices,2,Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis = (occupancy[cat2] - occupancy[cat1]) * log((*this)[cat1] / (*this)[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted)	{
                    total += 1.0;
                    occupancy.Swap(cat1,cat2);
                    permut.Swap(cat1,cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    double MoveAdjacentCompAlloc(int k0, OccupancySuffStat& occupancy, Permutation& permut)	{

        GibbsResample(occupancy);
        int nrep = (int) (k0 * kappa);
        
        double total = 0;

        for (int i=0; i<nrep; i++)	{
            int cat1 = (int)(Random::Uniform() * (GetSize()-2));  
            int cat2 = cat1 + 1;
            double logMetropolis = (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1-V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted)	{
                total += 1.0;
                SwapComponents(cat1,cat2);
                occupancy.Swap(cat1,cat2);
                permut.Swap(cat1,cat2);
            }
        }

        return total /= nrep;
    }

    void FromStreamSB(istream& is) {
        for (int k=0; k<GetSize(); k++)  {
            is >> (*this)[k];
            is >> V[k];
        }
    }

    void ToStreamSB(ostream& os) const {
        for (int k=0; k<GetSize(); k++)  {
            os << GetVal(k) << '\t';
            os << V[k] << '\t';
        }
    }

    unsigned int GetMPISizeSB() const  {
        return 2*GetSize();
    }

    //! get array from MPI buffer
    void MPIGetSB(const MPIBuffer& is)    {
        for (int k=0; k<GetSize(); k++) {
            is >> (*this)[k] >> V[k];
        }
    }

    //! write array into MPI buffer
    void MPIPutSB(MPIBuffer& os) const {
        for (int k=0; k<GetSize(); k++) {
            os << GetVal(k) << V[k];
        }
    }

    private:

    double kappa;
    vector<double> V;

};

#endif

