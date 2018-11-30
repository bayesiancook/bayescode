#include <sstream>
#include "mpi_components/PartitionedBufferManager.hpp"
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/reduce.hpp"
#include "operations/proxies.hpp"

using MPI::p;
using namespace std;

// clang-format off
struct MyStruct {
    int a;
    double b;
    double c;  // unused
    template <class T> void serialization_interface(T& x) { x.add(a, b); }
};
template<> struct has_custom_serialization<MyStruct> : public true_type {};
// clang-format on

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

// clang-format off
struct StructTheGreat {
    int a;
    double b;
    template <class T> void serialization_interface(T& x) { x.add(a, b); }
};
template<> struct has_custom_serialization<StructTheGreat> : public true_type {};
// clang-format on

void compute(int, char**) {
    if (!p->rank) {
        SendBuffer buf;
        int i[4] = {2, 3, 4, 5};
        double j[4] = {2, 3, 4, 5.2};
        buf.pack(i, 4);
        buf.pack(j, 4);

        ReceiveBuffer rcvbuf(buf.size());
        memcpy(rcvbuf.data(), buf.data(), buf.size());
        p->message("Unpacked ints: %d, %d, %d, %d", rcvbuf.unpack<int>(), rcvbuf.unpack<int>(),
            rcvbuf.unpack<int>(), rcvbuf.unpack<int>());
        p->message("Unpacked doubles: %.2f, %.2f, %.2f, %.2f", rcvbuf.unpack<double>(),
            rcvbuf.unpack<double>(), rcvbuf.unpack<double>(), rcvbuf.unpack<double>());

        ReceiveBuffer rcvbuf2(buf.size());
        memcpy(rcvbuf2.data(), buf.data(), buf.size());
        auto v1 = rcvbuf2.unpack_vector<int>(4);
        p->message("Unpacked ints: %d, %d, %d, %d", v1.at(0), v1.at(1), v1.at(2), v1.at(3));
        auto v2 = rcvbuf2.unpack_vector<double>(4);
        p->message(
            "Unpacked doubles: %.2f, %.2f, %.2f, %.2f", v2.at(0), v2.at(1), v2.at(2), v2.at(3));
    }
    if (!p->rank) {
        int i{1}, j{2}, k{3};
        double l{2.1}, m{3.2}, n{5.6};
        vector<int> v{5, 8, 9};
        vector<double> v2{2.35, 5.68};

        BufferManager b;
        b.add(i, j, k, l, m, n, v, v2);

        void* buf = b.send_buffer();
        size_t buf_size = b.buffer_size();
        p->message("Size of buffer is %d", buf_size);

        void* rcvbuf = b.receive_buffer();
        memcpy(rcvbuf, buf, buf_size);

        b.receive();
        p->message("%d, %d, %d, %.2f, %.2f, %.2f, {%d, %d, %d}, {%.2f, %.2f}", i, j, k, l, m, n,
            v.at(0), v.at(1), v.at(2), v2.at(0), v2.at(1));
    }
    if (!p->rank) {
        IndexSet is{"a", "b", "c", "d", "e"};
        Partition part(is, 3);
        PartitionedBufferManager bm(part);

        std::vector<double> v = {1.1, 2.2, 3.3, 4.4, 5.5};
        bm.add(v);

        std::vector<StructTheGreat> v2 = {{1, 0.1}, {2, 0.2}, {3, 0.3}, {4, 0.4}, {5, 0.5}};
        bm.add(v2);
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

    Group slave_operations{broadcast(m, {"a", "c"}), broadcast(m, {"i", "j"})};

    p->rank ? slave_operations.acquire() : slave_operations.release();
    p->rank ? master_operations.release() : master_operations.acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }