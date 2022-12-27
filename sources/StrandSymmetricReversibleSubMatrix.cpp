#include "StrandSymmetricReversibleSubMatrix.hpp"

void StrandSymmetricReversibleSubMatrix::ComputeStationary() const {
    mStationary[0] = (1.0 - gc) / 2;
    mStationary[1] = gc / 2;
    mStationary[2] = gc / 2;
    mStationary[3] = (1.0 - gc) / 2;
}

/* general strand symmetric - reversible
rho_AC rho_AG rho_AT rho_CG
*/

void StrandSymmetricReversibleSubMatrix::ComputeArray(int i) const {
    if (i == 0) {
        Q(0, 1) = rho[0] * gc / 2;
        Q(0, 2) = rho[1] * kappa * gc / 2;
        Q(0, 3) = rho[2] * (1.0 - gc) / 2;
        Q(0, 0) = -Q(0, 1) - Q(0, 2) - Q(0, 3);
    } else if (i == 1) {
        Q(1, 0) = rho[0] * (1.0 - gc) / 2;
        Q(1, 2) = rho[3] * gc / 2;
        Q(1, 3) = rho[1] * kappa * (1.0 - gc) / 2;
        Q(1, 1) = -Q(1, 0) - Q(1, 2) - Q(1, 3);
    } else if (i == 2) {
        Q(2, 0) = rho[1] * kappa * (1.0 - gc) / 2;
        Q(2, 1) = rho[3] * gc / 2;
        Q(2, 3) = rho[0] * (1.0 - gc) / 2;
        Q(2, 2) = -Q(2, 0) - Q(2, 1) - Q(2, 3);
    } else if (i == 3) {
        Q(3, 0) = rho[2] * (1 - 0 - gc) / 2;
        Q(3, 1) = rho[1] * kappa * gc / 2;
        Q(3, 2) = rho[0] * gc / 2;
        Q(3, 3) = -Q(3, 0) - Q(3, 1) - Q(3, 2);
    }
}
