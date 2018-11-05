#include <sstream>
#include "broadcast.hpp"
#include "gather.hpp"

using MPI::p;
using namespace std;

struct DummyModel {
    double a{-1}, b{-1}, c{-1};
    vector<double> v;
    Partition partition{{"d", "e", "f"}, static_cast<size_t>(p->size) - 1, 1};

    template <class C>
    void declare_model(C& ref) {
        ref.add("a", a);
        ref.add("b", b);
        ref.add("c", c);
        ref.add("v", v, partition);
    }

    void print() {
        stringstream ss;
        for (auto e : v) { ss << e << " "; }
        p->message("Model state is %f, %f, %f, {%s}", a, b, c, ss.str().c_str());
    }
};

void compute(int, char**) {
    DummyModel m;
    if (!p->rank) {  // master
        m.a = 2.2;
        m.b = 3.3;
        m.c = 4.4;
    } else {  // slave
        m.v = std::vector<double>(m.partition.my_partition_size(), p->rank + 1.1);
    }

    auto gatherer = gather_model(m, {"v"});
    auto bcaster = broadcast_model(m, {"a", "c"});

    p->rank ? bcaster->acquire() : bcaster->release();
    p->rank ? gatherer->release() : gatherer->acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }