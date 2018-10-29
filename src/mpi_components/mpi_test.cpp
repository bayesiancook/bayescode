#include "Process.hpp"

using namespace MPI;

void compute(int, char**) {
    p.message("Hello world!");
    MPI_Barrier(MPI_COMM_WORLD);
    p.message("Goodbye world!");
}

int main(int argc, char** argv) {
    mpi_run(argc, argv, compute);
}