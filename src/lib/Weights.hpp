#pragma once

#include "Array.hpp"
#include "OccupancySuffStat.hpp"
#include "Random.hpp"

/**
 * \brief A truncated stick-breaking weight vector
 */

class Weights : public SimpleArray<double> {
  public:
    //! constructor, parameterized by array size (truncation upper limit) and
    //! concentration parameter kappa
    Weights(int inncat, double &inp0) : SimpleArray<double>(inncat), p0(inp0) { Update(); }

    ~Weights() override = default;

    double GetWeight(int cat) {
        if (cat == 0) {
            return p0;
        } else {
            return (1 - p0) / (this->GetSize() - 1);
        }
    }

    void Update() {
        for (int cat{0}; cat < this->GetSize(); cat++) { (*this)[cat] = GetWeight(cat); }
    }

    void GibbsResample(OccupancySuffStat const &occupancy, double mean, double invconc) {
        double a = mean / invconc;
        double b = 1.0 / invconc - a;

        int n0 = occupancy.GetVal(0), n1 = 0;
        assert(occupancy.GetSize() == this->GetSize());
        for (int cat{1}; cat < this->GetSize(); cat++) { n1 += occupancy.GetVal(cat); }

        double x = Random::sGamma(a + n0);
        double y = Random::sGamma(b + n1);
        p0 = x / (x + y);
        Update();
    }

  private:
    double &p0;
};
