#pragma once

#include <assert.h>
#include "lib/Random.hpp"

class AcceptanceStats {
    double ntot{0};
    double nacc{0};

  public:
    void accept() {
        ntot++;
        nacc++;
    }
    void reject() { ntot++; }
    double ratio() const { return nacc / ntot; }
};

bool decide(double log_prob) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    return log(uniform(Random::global_gen)) < log_prob;
}