#include <mpi.h>
#include <cstdio>

int main() {
    // MPI init
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // data
    int a{-1};
    float b{-1};
    struct {
        int c{-1};
        float d{-1};
    } s;

    // init on sender
    if (rank == 0) {
        a = 3;
        b = 5.2;
        s.c = 9;
        s.d = 6.32;
    }

    // addresses
    MPI_Aint aa, ab, ac, ad;
    MPI_Get_address(&a, &aa);
    MPI_Get_address(&b, &ab);
    MPI_Get_address(&s.c, &ac);
    MPI_Get_address(&s.d, &ad);

    // mpi type data
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_INT, MPI_FLOAT, MPI_INT, MPI_FLOAT};
    MPI_Datatype mpi_car_type;
    MPI_Aint offsets[4] = {0, ab - aa, ac - aa, ad - aa};

    // create and commit type
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_car_type);
    MPI_Type_commit(&mpi_car_type);

    // actual comm + test
    printf("<%d> Data: %d, %f, %d, %f\n", rank, a, b, s.c, s.d);
    MPI_Bcast(&a, 1, mpi_car_type, 0, MPI_COMM_WORLD);
    printf("<%d> Data: %d, %f, %d, %f\n", rank, a, b, s.c, s.d);

    // MPI stuff
    MPI_Finalize();
}