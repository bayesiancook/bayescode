#include "StrandSymmetricIrreversibleSubMatrix.hpp"

void StrandSymmetricIrreversibleSubMatrix::ComputeStationary() const {
    mStationary[0] = 0.25;
    mStationary[1] = 0.25;
    mStationary[2] = 0.25;
    mStationary[3] = 0.25;
}

/* general strand symmetric - irreversible
AC - TG
AG - TC
AT - TA
CA - GT
CG - GC
CT - GA
*/

void StrandSymmetricIrreversibleSubMatrix::ComputeArray(int i) const {
    if (i == 0) {
        Q(0, 1) = rates[0];
        Q(0, 2) = rates[1];
        Q(0, 3) = rates[2];
        Q(0, 0) = -Q(0, 1) - Q(0, 2) - Q(0, 3);
    } else if (i == 1) {
        Q(1, 0) = rates[3];
        Q(1, 2) = rates[4];
        Q(1, 3) = rates[5];
        Q(1, 1) = -Q(1, 0) - Q(1, 2) - Q(1, 3);
    } else if (i == 2) {
        Q(2, 0) = rates[5];
        Q(2, 1) = rates[4];
        Q(2, 3) = rates[3];
        Q(2, 2) = -Q(2, 0) - Q(2, 1) - Q(2, 3);
    } else if (i == 3) {
        Q(3, 0) = rates[2];
        Q(3, 1) = rates[1];
        Q(3, 2) = rates[0];
        Q(3, 3) = -Q(3, 0) - Q(3, 1) - Q(3, 2);
    }
}

/*
void StrandSymmetricIrreversibleSubMatrix::ComputeArray(int i) const {
    if (i == 0) {
        Q(0, 1) = rates(0);
        Q(0, 2) = rates(1);
        Q(0, 3) = rates(2);
        Q(0, 0) = -Q(0, 1) - Q(0, 2) - Q(0, 3);
    } else if (i == 1) {
        Q(1, 0) = rates(3);
        Q(1, 2) = rates(4);
        Q(1, 3) = rates(5);
        Q(1, 1) = -Q(1, 0) - Q(1, 2) - Q(1, 3);
    } else if (i == 2) {
        Q(2, 0) = rates(5);
        Q(2, 1) = rates(4);
        Q(2, 3) = rates(3);
        Q(2, 2) = -Q(2, 0) - Q(2, 1) - Q(2, 3);
    } else if (i == 3) {
        Q(3, 0) = rates(2);
        Q(3, 1) = rates(1);
        Q(3, 2) = rates(0);
        Q(3, 3) = -Q(3, 0) - Q(3, 1) - Q(3, 2);
    }
}
*/
