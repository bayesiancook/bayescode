#include <sstream>
#include "doctest.h"
#include "mpi_components/partition.hpp"
#include "mpi_components/utils.hpp"

using namespace std;

TEST_CASE("Full Partition test") {
    IndexSet s{"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"};
    IndexSet ss1{"a", "b", "c", "d"};
    IndexSet ss2{"e", "f", "g", "h", "i"};
    IndexSet ss3{"j", "k", "l", "m", "n"};

    // partition as if there was 1 master and 3 slaves
    Process process0(0, 4);
    Process process1(1, 4);
    Process process2(2, 4);
    Process process3(3, 4);
    Partition p0(s, 3, 1, process0);
    Partition p1(s, 3, 1, process1);
    Partition p2(s, 3, 1, process2);
    Partition p3(s, 3, 1, process3);

    // sizes
    CHECK(p0.size_all() == 14);
    CHECK(p1.size_all() == 14);
    CHECK(p2.size_all() == 14);
    CHECK(p3.size_all() == 14);
    CHECK(p0.size() == 3);
    CHECK(p1.size() == 3);
    CHECK(p2.size() == 3);
    CHECK(p3.size() == 3);
    CHECK(p1.my_partition_size() == 4);
    CHECK(p2.my_partition_size() == 5);
    CHECK(p3.my_partition_size() == 5);
    CHECK(p0.max_partition_size() == 5);
    CHECK(p1.max_partition_size() == 5);
    CHECK(p2.max_partition_size() == 5);
    CHECK(p3.max_partition_size() == 5);
    CHECK(p0.partition_size(1) == 4);
    CHECK(p0.partition_size(2) == 5);
    CHECK(p0.partition_size(3) == 5);
    CHECK(p2.partition_size(1) == 4);
    CHECK(p2.partition_size(2) == 5);
    CHECK(p2.partition_size(3) == 5);

    // partitions
    CHECK(p1.my_partition() == ss1);
    CHECK(p2.my_partition() == ss2);
    CHECK(p3.my_partition() == ss3);
    CHECK(p0.get_all() == s);
    CHECK(p2.get_all() == s);
    CHECK(p0.get_partition(1) == ss1);
    CHECK(p0.get_partition(3) == ss3);
    CHECK(p2.get_partition(1) == ss1);
    CHECK(p2.get_partition(3) == ss3);

    // owner
    CHECK(p0.owner("a") == 1);
    CHECK(p0.owner("e") == 2);
    CHECK(p0.owner("k") == 3);
    CHECK(p2.owner("a") == 1);
    CHECK(p2.owner("e") == 2);
    CHECK(p2.owner("k") == 3);
}

TEST_CASE("Partition: obtain local process") {
    IndexSet s{"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"};
    IndexSet ss2{"e", "f", "g", "h", "i"};

    // partition as if there was 1 master and 3 slaves
    MPI::p = std::unique_ptr<Process>(new Process(2, 4));
    Partition p(s, 3, 1);
    CHECK(p.my_partition_size() == 5);
    CHECK(p.my_partition() == ss2);
}

TEST_CASE("Process::message") {
    stringstream ss;
    Process p(1, 2);
    p.message(ss, "Hello %d - %s", 3, "aa");

    CHECK(ss.str() ==
          "\e[0m\e[1m[\e[0m\e[32m1\e[0m\e[1m/\e[0m\e[32m2\e[0m\e[1m] \e[0mHello 3 - aa\n");
}