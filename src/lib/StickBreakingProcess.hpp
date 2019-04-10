#ifndef SBP_H
#define SBP_H

#include "Array.hpp"
#include "OccupancySuffStat.hpp"
#include "Permutation.hpp"
#include "global/Random.hpp"

/**
 * \brief A truncated stick-breaking weight vector
 */

class StickBreakingProcess : public SimpleArray<double> {
  public:
    //! constructor, parameterized by array size (truncation upper limit) and
    //! concentration parameter kappa
    StickBreakingProcess(int inncat, const double &inkappa)
        : SimpleArray<double>(inncat), V(inncat), kappa(inkappa) {
        Sample();
    }

    ~StickBreakingProcess() {}

    const double &GetKappa() const { return kappa; }

    //! get underlying array of Beta variates
    std::vector<double> &GetBetaVariates() { return V; }
    const std::vector<double> &GetBetaVariates() const { return V; }

    //! swap the two components
    void SwapComponents(int cat1, int cat2) {
        Swap(cat1, cat2);
        double tmp = V[cat1];
        V[cat1] = V[cat2];
        V[cat2] = tmp;
    }

    //! Sample (from prior distribution)
    void Sample() {
        double cumulProduct = 1.0;
        double totweight = 0;
        for (int k = 0; k < GetSize(); k++) {
            double x = Random::sGamma(1.0);
            double y = Random::sGamma(kappa);
            double v = x / (x + y);
            V[k] = v;
            if (k == GetSize() - 1) {
                V[k] = 1;
                v = 1;
            }
            (*this)[k] = v * cumulProduct;
            cumulProduct *= (1 - v);
            totweight += (*this)[k];
        }
    }

    //! Sample from posterior distribution (given occupany sufficient statistics)
    void GibbsResample(const OccupancySuffStat &occupancy) {
        int remainingOcc = 0;
        for (int i = 0; i < occupancy.GetSize(); i++) { remainingOcc += occupancy.GetVal(i); }

        double cumulProduct = 1.0;
        double totweight = 0;
        for (int k = 0; k < GetSize(); k++) {
            remainingOcc -= occupancy.GetVal(k);
            double x = Random::sGamma(1 + occupancy.GetVal(k));
            double y = Random::sGamma(kappa + remainingOcc);
            double v = x / (x + y);
            if (!v) { v = 1e-50; }
            V[k] = v;
            if (k == GetSize() - 1) {
                V[k] = 1;
                v = 1;
            }
            (*this)[k] = v * cumulProduct;
            cumulProduct *= (1 - v);
            totweight += (*this)[k];
        }
    }

    //! get log prior (under current kappa value)
    double GetLogProb() const {
        double total = 0;
        for (int k = 0; k < GetSize() - 1; k++) {
            total += Random::logBetaDensity(V[k], 1.0, kappa);
        }
        return total;
    }

    //! get log prior (under specified kappa value)
    double GetLogProb(double kappa) const {
        double total = 0;
        for (int k = 0; k < GetSize() - 1; k++) {
            total += Random::logBetaDensity(V[k], 1.0, kappa);
        }
        return total;
    }

    //! get marginal log prior (i.e. integrated over weights, conditional on
    //! specified occupancy sufficient statistics)
    double GetMarginalLogProb(const OccupancySuffStat &occupancy) const {
        int remainingOcc = 0;
        for (int i = 0; i < occupancy.GetSize(); i++) { remainingOcc += occupancy.GetVal(i); }

        double total = 0;
        for (int k = 0; k < GetSize(); k++) {
            if (remainingOcc) {
                remainingOcc -= occupancy.GetVal(k);
                total += log(kappa) + Random::logGamma(1 + occupancy.GetVal(k)) +
                         Random::logGamma(kappa + remainingOcc) -
                         Random::logGamma(1 + kappa + occupancy.GetVal(k) + remainingOcc);
            }
        }
        if (remainingOcc) {
            std::cerr << "error in allocation count\n";
            exit(1);
        }
        return total;
    }

    //! perform nrep cycles of Label-switching MCMC updates (such as defined in
    //! Papaspiliopoulos and Roberts, 2008, Biometrika, 95:169)
    void LabelSwitchingMove(int nrep, OccupancySuffStat &occupancy, Permutation &permut) {
        MoveOccupiedCompAlloc(nrep, occupancy, permut);
        MoveAdjacentCompAlloc(nrep, occupancy, permut);
    }

    //! label-switching MCMC update: swapping two randomly chosen occupied
    //! components
    double MoveOccupiedCompAlloc(int k0, OccupancySuffStat &occupancy, Permutation &permut) {
        int nrep = (int)(k0 * kappa);
        GibbsResample(occupancy);
        double total = 0.0;
        int Nocc = occupancy.GetNcluster();
        if (Nocc != 1) {
            for (int i = 0; i < nrep; i++) {
                int occupiedComponentIndices[Nocc];
                int j = 0;
                for (int k = 0; k < GetSize(); k++) {
                    if (occupancy[k] != 0) {
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc) {
                    std::cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices, 2, Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis =
                    (occupancy[cat2] - occupancy[cat1]) * log((*this)[cat1] / (*this)[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted) {
                    total += 1.0;
                    occupancy.Swap(cat1, cat2);
                    permut.Swap(cat1, cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    //! label-switching MCMC udpate: swapping two successive components
    double MoveAdjacentCompAlloc(int k0, OccupancySuffStat &occupancy, Permutation &permut) {
        GibbsResample(occupancy);
        int nrep = (int)(k0 * kappa);

        double total = 0;

        for (int i = 0; i < nrep; i++) {
            int cat1 = (int)(Random::Uniform() * (GetSize() - 2));
            int cat2 = cat1 + 1;
            double logMetropolis =
                (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1 - V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted) {
                total += 1.0;
                SwapComponents(cat1, cat2);
                occupancy.Swap(cat1, cat2);
                permut.Swap(cat1, cat2);
            }
        }

        return total /= nrep;
    }

    template <class Info>
    void declare_interface(Info info) {
        declare(info, "array", array);
        declare(info, "betavariates", GetBetaVariates());
    }

  private:
    std::vector<double> V;
    const double &kappa;
};

#endif
