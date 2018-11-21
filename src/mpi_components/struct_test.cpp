#include <mpi.h>
#include <cstdio>
#include "utils.hpp"

MPI_Datatype MPI_MYSTRUCT;

int main() {
    // MPI init
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // data
    struct Mystruct {
        int a{-1};
        float b{-1};
    };
    Mystruct s;

    // init on sender
    if (rank == 0) {
        s.a = 9;
        s.b = 6.32;
    }

    STRUCT_DECL(Mystruct)
    ATTRIBUTE(a)
    ATTRIBUTE(b)
    STRUCT_DECL_COMMIT(MPI_MYSTRUCT)

    // actual comm + test
    printf("<%d> Data: %d, %f\n", rank, s.a, s.b);
    MPI_Bcast(&s, 1, MPI_MYSTRUCT, 0, MPI_COMM_WORLD);
    printf("<%d> Data: %d, %f\n", rank, s.a, s.b);

    // MPI stuff
    MPI_Finalize();
}