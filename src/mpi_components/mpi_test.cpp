#include "Broadcaster.hpp"

using MPI::p;
using namespace std;

struct DummyModel {
    double a, b, c;
    vector<double> v;

    template <class C>
    void declare_model(C& ref) {
        ref.add("a", a);
        ref.add("b", b);
        ref.add("c", c);
        ref.add("v", v);
    }

    void print() {
        p->message("Model state is %f, %f, %f, {%f, %f, %f}", a, b, c, v[0], v[1], v[2]);
    }
};

void compute(int, char**) {
    DummyModel m = {-1, -1, -1, {-1, -1, -1}};
    if (!p->rank) { m = {12.2, 15.5, 16.8, {7.7, 8.8, 9.9}}; }

    auto bcaster = broadcast_model(m, {"a", "c", "v"});

    if (p->rank) { bcaster->acquire(); }
    bcaster->release();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }