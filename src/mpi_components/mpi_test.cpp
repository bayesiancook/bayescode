#include "Broadcaster.hpp"

using MPI::p;

void compute(int, char**) {
    p->message("Hello world!");
    MPI_Barrier(MPI_COMM_WORLD);
    p->message("Goodbye world!");


    if (!p->rank) {
        BroadcasterMaster<double> bcast;
        double a{3.1}, b{2.3}, c{4.5};
        bcast.add(a);
        bcast.add(b);
        bcast.add(c);
        bcast.release();
    } else {
        BroadcasterSlave<double> bcast;
        double a{-1}, b{-1}, c{-1};
        bcast.add(a);
        bcast.add(b);
        bcast.add(c);
        bcast.acquire();
        p->message("Got values %f, %f and %f", a, b, c);
    }
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }