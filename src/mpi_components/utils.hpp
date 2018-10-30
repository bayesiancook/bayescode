#pragma once

#include <mpi.h>

template <typename T>
MPI_Datatype get_datatype() {
    if (std::is_same<T, double>()) {
        return MPI_DOUBLE;
    } else if (std::is_same<T, int>()) {
        return MPI_INT;
    } else {
        fprintf(
            stderr, "Error in get_datatype: template parameter T is not a supported datatype!\n");
        exit(1);
    }
}