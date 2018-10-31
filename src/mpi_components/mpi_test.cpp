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

    void print() { p->message("Model state is %f, %f, %f", a, b, c); }
};

void compute(int, char**) {
    DummyModel m = {-1, -1, -1};
    if (!p->rank) { m = {12.2, 15.5, 16.8}; }

    auto bcaster = broadcast_model(m, {"a", "c"});

    if (p->rank) { bcaster->acquire(); }
    bcaster->release();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }