#include <sstream>
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/reduce.hpp"
#include "mpi_components/serialization.hpp"
#include "operations/proxies.hpp"

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
        p->message("Model state is %.2f, %.2f, %.2f, {%s}, %.2f, %.2f, (%d, %.2f), (%d, %.2f)", a,
            b, c, ss.str().c_str(), g, h, i.a, i.b, j.a, j.b);
    }
};

void compute(int, char**) {
    STRUCT_DECL(MyStruct)
    ATTRIBUTE(a)
    ATTRIBUTE(b)
    STRUCT_COMMIT(MPI_MYSTRUCT)

    if (!p->rank) {
        SendBuffer buf;
        int i[4] = {2, 3, 4, 5};
        double j[4] = {2, 3, 4, 5.2};
        buf.pack(i, 4);
        buf.pack(j, 4);

        ReceiveBuffer rcvbuf(buf.data(), buf.size());
        p->message("Unpacked ints: %d, %d, %d, %d", rcvbuf.unpack<int>(), rcvbuf.unpack<int>(),
            rcvbuf.unpack<int>(), rcvbuf.unpack<int>());
        p->message("Unpacked doubles: %.2f, %.2f, %.2f, %.2f", rcvbuf.unpack<double>(),
            rcvbuf.unpack<double>(), rcvbuf.unpack<double>(), rcvbuf.unpack<double>());
    }

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

    Group master_operations{reduce<double>(m, {"g", "h"}), gather<double>(m, {"v"})};

    Group slave_operations{broadcast<double>(m, {"a", "c"}), broadcast<MyStruct>(m, {"i", "j"})};

    p->rank ? slave_operations.acquire() : slave_operations.release();
    p->rank ? master_operations.release() : master_operations.acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }