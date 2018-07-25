#include "T92SubMatrix.hpp"

void T92SubMatrix::ComputeStationary() const {
    mStationary[0] = (1.0 - gc) / 2;
    mStationary[1] = gc / 2;
    mStationary[2] = gc / 2;
    mStationary[3] = (1.0 - gc) / 2;
}

void T92SubMatrix::ComputeArray(int i) const {
    if (i == 0) {
        Q(0, 1) = gc / 2;
        Q(0, 2) = kappa * gc / 2;
        Q(0, 3) = (1.0 - gc) / 2;
        Q(0, 0) = -Q(0, 1) - Q(0, 2) - Q(0, 3);
    } else if (i == 1) {
        Q(1, 0) = (1.0 - gc) / 2;
        Q(1, 2) = gc / 2;
        Q(1, 3) = kappa * (1.0 - gc) / 2;
        Q(1, 1) = -Q(1, 0) - Q(1, 2) - Q(1, 3);
    } else if (i == 2) {
        Q(2, 0) = kappa * (1.0 - gc) / 2;
        Q(2, 1) = gc / 2;
        Q(2, 3) = (1.0 - gc) / 2;
        Q(2, 2) = -Q(2, 0) - Q(2, 1) - Q(2, 3);
    } else if (i == 3) {
        Q(3, 0) = (1 - 0 - gc) / 2;
        Q(3, 1) = kappa * gc / 2;
        Q(3, 2) = gc / 2;
        Q(3, 3) = -Q(3, 0) - Q(3, 1) - Q(3, 2);
    }
}
