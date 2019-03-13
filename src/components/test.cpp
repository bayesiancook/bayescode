#include "doctest.h"

#include "BaseArgParse.hpp"
#include "Tracer.hpp"

using namespace std;

struct TracerDummyModel {
    double a, b;
    vector<double> v;

    template <class Info>
    void declare_interface(Info info) {
        declare(info, "a", a);
        declare(info, "v", v);
        declare(info, "b", b);
    }
};

TEST_CASE("Tracer reading test") {
    stringstream ss("1.1 2.2 3.3 4.4 5.5");
    TracerDummyModel m = {-1, -1, {-1, -1, -1}};
    Tracer t(m, processing::True());

    t.read_line(ss);

    CHECK(m.a == 1.1);
    CHECK(m.b == 5.5);
    vector<double> expected{2.2, 3.3, 4.4};
    CHECK(m.v == expected);
}

TEST_CASE("Tracer writing test") {
    stringstream ss;
    TracerDummyModel m = {1.1, 2.2, {3.3, 4.4, 5.5}};
    Tracer t(m, processing::True());

    t.write_header(ss);
    t.write_line(ss);
    auto old_v_addr = &(m.v[0]);
    m.a = 1.11;
    m.b = 2.21;
    m.v[0] = 3.31;
    m.v[1] = 4.41;
    m.v[2] = 5.51;
    CHECK(old_v_addr == &(m.v[0]));
    CHECK(m.v[0] == 3.31);
    CHECK(m.v[1] == 4.41);
    CHECK(m.v[2] == 5.51);
    t.write_line(ss);

    CHECK(ss.str() ==
          "a\tv[0]\tv[1]\tv[2]\tb\n1.1\t3.3\t4.4\t5.5\t2.2\n1.11\t3.31\t4.41\t5.51\t2.21");
}

struct MyTracerData {
    int x, y;
    template <class Info>
    void declare_interface(Info info) {
        declare(info, "x", x);
        declare(info, "y", y);
    }
};

struct MyTracerStruct {
    int a{17};
    vector<int> b{2, 3, 4, 5};
    vector<MyTracerData> c{{1, 2}, {2, 3}};
    template <class Info>
    void declare_interface(Info info) {
        declare(info, "a", a);
        declare(info, "b", b);
        declare(info, "c", c);
    }
};

TEST_CASE("Tracer unrolling test") {
    MyTracerStruct s;
    Tracer t(s, processing::True());
    stringstream ss;

    t.write_header(ss);
    t.write_line(ss);
    CHECK(ss.str() ==
          "a	b[0]	b[1]	b[2]	b[3]	c_0_x	c_0_y	c_1_x	c_1_y\n17	2	"
          "3	4	5	1	2	2	3");
}

struct MyArgs : public BaseArgParse {
    MyArgs(ChainCmdLine& cmd) : BaseArgParse(cmd) {}
    ValueArg<std::string> treefile{"t", "tree", "", true, "", "string", cmd};
    ValueArg<int> until{"u", "until", "", false, -1, "int", cmd};
    SwitchArg force{"f", "force", "", cmd};
};

TEST_CASE("Arg parse test") {
    vector<string> argv_str = {"test_bin", "-t", "tree.tree", "-u", "19", "-f", "tmp"};
    int argc = argv_str.size();
    char* argv[argc];
    for (int i = 0; i < argc; i++) {
        argv[i] = static_cast<char*>(malloc(sizeof(char) * argv_str.at(i).size() + 1));
        strcpy(argv[i], argv_str.at(i).c_str());
    }

    ChainCmdLine cmd{argc, argv, "test_bin", ' ', "0.1"};
    MyArgs args(cmd);
    cmd.parse();

    CHECK(args.until.getValue() == 19);
    CHECK(args.force.getValue() == true);
    CHECK(args.treefile.getValue() == "tree.tree");
    CHECK(cmd.chain_name() == "tmp");
    // CHECK(cmd.checkpoint_file() == "tmp.param");

    for (int i = 0; i < argc; i++) { free(argv[i]); }
}