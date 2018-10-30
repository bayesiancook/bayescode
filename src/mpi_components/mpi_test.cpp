#include "Broadcaster.hpp"

using MPI::p;

struct DummyModel {
    double a, b, c;

    template <class C>
    void declare_model(C& ref) {
        ref.add("a", a);
        ref.add("b", b);
        ref.add("c", c);
    }
};

void compute(int, char**) {
    if (!p->rank) {
        DummyModel m = {2.3, 2.4, 5.2};
        BroadcasterMaster<double> bcast({"a", "c"});
        bcast.register_from_method(&m, &DummyModel::declare_model);
        bcast.release();
    } else {
        DummyModel m = {-1, -1, -1};
        BroadcasterSlave<double> bcast({"a", "c"});
        bcast.register_from_method(&m, &DummyModel::declare_model);
        bcast.acquire();
        p->message("Got values %f, %f and %f", m.a, m.b, m.c);
    }
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }