#pragma once

#include "Array.hpp"
#include "IncompleteGamma.hpp"

class DGamma : public SimpleArray<double> {
  public:
    DGamma(int Ncat, double const& mean, double const& invshape, double delta_shift)
        : SimpleArray<double>(Ncat), mean(mean), invshape(invshape) {
        // The first category is user defined (and not updated)
        (*this)[0] = delta_shift;
        Update();
    }

    ~DGamma() override = default;

    std::vector<double> GetDiscreteArray(int Ncat) {
        double alpha = 1.0 / invshape;
        std::vector<double> x(Ncat, 0.0);
        if (Ncat == 1) {
            x[0] = mean;
        } else {
            std::vector<double> y(Ncat, 0.0);
            double lg = Random::logGamma(alpha + 1.0);
            for (int i = 0; i < Ncat; i++) { x[i] = PointGamma((i + 1.0) / Ncat, alpha, alpha); }
            for (int i = 0; i < Ncat - 1; i++) {
                y[i] = IncompleteGamma(alpha * x[i], alpha + 1, lg);
            }
            y[Ncat - 1] = 1.0;

            x[0] = mean * Ncat * y[0];
            for (int i = 1; i < Ncat; i++) { x[i] = mean * Ncat * (y[i] - y[i - 1]); }
        }
        return x;
    }


    void Update() {
        auto x = GetDiscreteArray(this->GetSize() - 1);
        for (int i = 1; i < this->GetSize(); i++) { (*this)[i] = x[i - 1]; }
    }

    double const& mean;
    double const& invshape;
};
