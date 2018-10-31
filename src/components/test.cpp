#include "doctest.h"

#include "BaseArgParse.hpp"
#include "Tracer.hpp"

using namespace std;

struct TracerDummyModel {
    double a, b;
    vector<double> v;

    template <class T>
    void my_declare_model(T& m) {
        m.add("a", a);
        m.add("v", v);
        m.add("b", b);
    }
};

TEST_CASE("Tracer reading test") {
    stringstream ss("1.1 2.2 3.3 4.4 5.5");
    TracerDummyModel m = {-1, -1, {-1, -1, -1}};
    Tracer t(m, &TracerDummyModel::my_declare_model);

    t.read_line(ss);

    CHECK(m.a == 1.1);
    CHECK(m.b == 5.5);
    vector<double> expected{2.2, 3.3, 4.4};
    CHECK(m.v == expected);
}

TEST_CASE("Tracer writing test") {
    stringstream ss;
    TracerDummyModel m = {1.1, 2.2, {3.3, 4.4, 5.5}};
    Tracer t(m, &TracerDummyModel::my_declare_model);

    t.write_header(ss);
    t.write_line(ss);
    m = {1.11, 2.21, {3.31, 4.41, 5.51}};
    t.write_line(ss);

    CHECK(ss.str() ==
          "a\tv[0]\tv[1]\tv[2]\tb\n1.1\t3.3\t4.4\t5.5\t2.2\n1.11\t3.31\t4.41\t5.51\t2.21");
}

struct MyArgs : public BaseArgParse {};

TEST_CASE("Arg parse test") {}