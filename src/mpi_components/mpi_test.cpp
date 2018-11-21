#include <sstream>
#include "broadcast.hpp"
#include "gather.hpp"
#include "reduce.hpp"

using MPI::p;
using namespace std;

struct MyStruct {
    int a;
    float b;
};

MPI_Datatype MPI_MYSTRUCT;
template <>
struct GetMPIDatatype<MyStruct> {
    MPI_Datatype operator()() { return MPI_MYSTRUCT; }
};

struct DummyModel {
    double a{-1}, b{-1}, c{-1};
    vector<double> v;
    Partition partition{{"d", "e", "f"}, static_cast<size_t>(p->size) - 1, 1};
    double g{-1}, h{-1};
    MyStruct i{-1, -1};

    template <class C>
    void declare_model(C& ref) {
        ref.add("a", a);
        ref.add("b", b);
        ref.add("c", c);
        ref.add("v", v, partition);
        ref.add("g", g);
        ref.add("h", h);
        ref.add("i", i);
    }

    void print() {
        stringstream ss;
        for (auto e : v) { ss << e << " "; }
        p->message("Model state is %f, %f, %f, {%s}, %f, %f, (%d, %f)", a, b, c, ss.str().c_str(),
            g, h, i.a, i.b);
    }
};

void compute(int, char**) {
    STRUCT_DECL(MyStruct)
    ATTRIBUTE(a)
    ATTRIBUTE(b)
    STRUCT_DECL_COMMIT(MPI_MYSTRUCT)

    DummyModel m;
    if (!p->rank) {  // master
        m.a = 2.2;
        m.b = 3.3;
        m.c = 4.4;
        m.i = {2, 3.2};
    } else {  // slave
        m.v = std::vector<double>(m.partition.my_partition_size(), p->rank + 1.1);
        m.g = p->rank + 0.4;
        m.h = p->rank - 0.4;
    }

    auto reducer = reduce_model(m, {"g", "h"});
    auto gatherer = gather_model(m, {"v"});
    auto bcaster = broadcast_model(m, {"a", "c"});
    auto struct_bcaster = broadcast_model<DummyModel, MyStruct>(m, {"i"});

    p->rank ? bcaster->acquire() : bcaster->release();
    p->rank ? struct_bcaster->acquire() : struct_bcaster->release();
    p->rank ? gatherer->release() : gatherer->acquire();
    p->rank ? reducer->release() : reducer->acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }