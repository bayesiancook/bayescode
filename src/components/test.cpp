#include "doctest.h"

#include "BaseArgParse.hpp"
#include "Tracer.hpp"
#include "monitoring.hpp"

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

class MyMonitor : public AbstractMonitor {
    int count;

  public:
    MyMonitor(int count) : count(count) {}

    void print(std::ostream& os) const final { os << "MyMonitor(" << count << ")"; }

    void update(int n) { count += n; }
};

TEST_CASE("Monitoring base classes") {
    auto f = [](int m) { return 3 * m; };

    MonitorManager m;
    m.new_monitor<MyMonitor>("m1", 3);
    m.new_monitor<MyMonitor>("m2", 5);

    std::stringstream ss;
    m.print(ss);
    CHECK(ss.str() == "m1: MyMonitor(3)\nm2: MyMonitor(5)\n");
    ss.str("");

    m.run_and_monitor<MyMonitor>("m1", f, 4);
    m.run_and_monitor<MyMonitor>("m2", f, 2);
    m.print(ss);
    CHECK(ss.str() == "m1: MyMonitor(15)\nm2: MyMonitor(11)\n");
}

TEST_CASE("Global monitor") {
    auto f = [](int m) { return 3 * m; };

    {
        gm->new_monitor<MyMonitor>("m1", 3);
        gm->new_monitor<MyMonitor>("m2", 5);
    }
    {
        gm->run_and_monitor<MyMonitor>("m1", f, 4);
        gm->run_and_monitor<MyMonitor>("m2", f, 2);
    }

    std::stringstream ss;
    gm->print(ss);
    CHECK(ss.str() == "m1: MyMonitor(15)\nm2: MyMonitor(11)\n");
}