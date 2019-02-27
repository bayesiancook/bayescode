#pragma once

#include "Array.hpp"
#include "OccupancySuffStat.hpp"
#include "Permutation.hpp"
#include "global/Random.hpp"

/**
 * \brief A truncated stick-breaking weight vector
 */

class Dirichlet : public SimpleArray<double> {
  public:
    //! constructor, parameterized by array size (truncation upper limit) and
    //! concentration parameter kappa
    Dirichlet(int inncat, const double &inkappa) : SimpleArray<double>(inncat), kappa(inkappa) {
        Sample();
    }

    ~Dirichlet() = default;

    const double &GetKappa() const { return kappa; }

    //! Sample (from prior distribution)
    void Sample() {
        double totweight = 0;
        for (int k = 0; k < GetSize(); k++) {
            (*this)[k] = Random::sGamma(kappa);
            totweight += (*this)[k];
        }
        for (int k = 0; k < GetSize(); k++) { (*this)[k] /= totweight; }
    }

    //! Sample from posterior distribution (given occupany sufficient statistics)
    void GibbsResample(const OccupancySuffStat &occupancy) {
        double totweight = 0;
        for (int k = 0; k < GetSize(); k++) {
            (*this)[k] = Random::sGamma(kappa + occupancy.GetVal(k));
            totweight += (*this)[k];
        }
        for (int k = 0; k < GetSize(); k++) { (*this)[k] /= totweight; }
    }

    //! get log prior (under current kappa value)
    double GetLogProb() const { return GetLogProb(kappa); }

    //! get log prior (under specified kappa value)
    double GetLogProb(double kappa) const {
        double tot = 0;
        for (int k = 0; k < GetSize(); k++) { tot += (kappa - 1) * log(GetVal(k)); }
        tot += -Random::logGamma(kappa) * GetSize() + Random::logGamma(kappa * GetSize());
        return tot;
    }

    //! get marginal log prior (i.e. integrated over weights, conditional on
    //! specified occupancy sufficient statistics)
    double GetMarginalLogProb(const OccupancySuffStat &occupancy) const {
        double tot = 0;
        int tot_occ = 0;
        for (int k = 0; k < GetSize(); k++) {
            tot_occ += occupancy.GetVal(k);
            tot += Random::logGamma(kappa + occupancy.GetVal(k));
        }
        tot += -Random::logGamma(kappa) * GetSize() + Random::logGamma(kappa * GetSize());
        tot -= Random::logGamma(tot_occ + kappa * GetSize());
        return tot;
    }

  private:
    const double &kappa;
};
