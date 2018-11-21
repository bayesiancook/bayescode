#include <sstream>
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/reduce.hpp"

using MPI::p;
using namespace std;

struct MyStruct {
    int a;
    float b;
    double c;  // unused
};

STRUCT_GLOBAL_DECL(MyStruct, MPI_MYSTRUCT)

struct DummyModel {
    double a{-1}, b{-1}, c{-1};
    vector<double> v;
    Partition partition{{"d", "e", "f"}, static_cast<size_t>(p->size) - 1, 1};
    MyStruct j{-1, -1, -1};
    double g{-1}, h{-1};
    MyStruct i{-1, -1, -1};

    template <class C>
    void declare_model(C& ref) {
        ref.add("a", a);
        ref.add("b", b);
        ref.add("c", c);
        ref.add("v", v, partition);
        ref.add("g", g);
        ref.add("h", h);
        ref.add("i", i);
        ref.add("j", j);
    }

    void print() {
        stringstream ss;
        for (auto e : v) { ss << e << " "; }
        p->message("Model state is %f, %f, %f, {%s}, %f, %f, (%d, %f), (%d, %f)", a, b, c,
            ss.str().c_str(), g, h, i.a, i.b, j.a, j.b);
    }
};

void compute(int, char**) {
    STRUCT_DECL(MyStruct)
    ATTRIBUTE(a)
    ATTRIBUTE(b)
    STRUCT_COMMIT(MPI_MYSTRUCT)

    DummyModel m;
    if (!p->rank) {  // master
        m.a = 2.2;
        m.b = 3.3;
        m.c = 4.4;
        m.i = {2, 3.2, -1};
        m.j = {7, 2.21, -1};
    } else {  // slave
        m.v = std::vector<double>(m.partition.my_partition_size(), p->rank + 1.1);
        m.g = p->rank + 0.4;
        m.h = p->rank - 0.4;
    }

    auto reducer = reduce_model(m, {"g", "h"});
    auto gatherer = gather_model(m, {"v"});
    auto bcaster = broadcast_model(m, {"a", "c"});
    auto struct_bcaster = broadcast_model<DummyModel, MyStruct>(m, {"i", "j"});

    p->rank ? bcaster->acquire() : bcaster->release();
    p->rank ? struct_bcaster->acquire() : struct_bcaster->release();
    p->rank ? gatherer->release() : gatherer->acquire();
    p->rank ? reducer->release() : reducer->acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }