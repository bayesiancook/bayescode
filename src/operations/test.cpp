#include "doctest.h"

#include <sstream>
#include "proxies.hpp"

using namespace std;

class MyProxy : public Proxy {
    int& target;

  public:
    MyProxy(int& target) : target(target) {}
    void acquire() final { target += 1; }
    void release() final { target *= 2; }
};

TEST_CASE("ForAll test") {
    int a = 1;
    int b = 1;
    MyProxy p1(a);
    MyProxy p2(b);
    MyProxy p3(a);
    CHECK(a == 1);
    CHECK(b == 1);

    ForAll group(&p1, &p2, &p3);

    group.acquire();
    CHECK(a == 3);
    CHECK(b == 2);

    group.release();
    CHECK(a == 12);
    CHECK(b == 4);

    group.add(&p3);  // re-adding p3 because why not

    group.acquire();
    CHECK(a == 15);
    CHECK(b == 5);

    group.release();
    CHECK(a == 120);
    CHECK(b == 10);
}

TEST_CASE("Group test") {
    int a = 1;
    int b = 1;
    std::unique_ptr<Proxy> p1(dynamic_cast<Proxy*>(new MyProxy(b)));
    // clang-format off
    Group group(
        std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new MyProxy(a))),
        std::move(p1),
        std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new MyProxy(b)))
    );
    // clang-format on
    CHECK(p1.get() == nullptr);

    group.acquire();  // 2 b and 1 a
    group.release();
    CHECK(a == 4);
    CHECK(b == 12);

    p1.reset(dynamic_cast<Proxy*>(new MyProxy(a)));
    group.add(std::move(p1));
    CHECK(p1.get() == nullptr);
    group.add(std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new MyProxy(b))));

    group.acquire();  // 3 b and 2 a
    group.release();
    CHECK(a == 24);
    CHECK(b == 120);
}

TEST_CASE("make_* functions") {
    int a = 1;
    int b = 1;
    std::unique_ptr<Proxy> p1(dynamic_cast<Proxy*>(new MyProxy(b)));
    MyProxy p3(a);
    MyProxy p4(b);

    // clang-format off
    auto group = make_group(
        std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new MyProxy(b))),
        std::move(p1)
    );
    auto forall = make_forall(&p3, &p4, group.get());
    // clang-format on
    CHECK(p1.get() == nullptr);

    forall->acquire();  // 3 b and 1 a
    forall->release();
    CHECK(a == 4);
    CHECK(b == 32);
}

TEST_CASE("Operation") {
    std::stringstream ss;
    Operation op([&ss]() { ss << "a = "; }, [&ss]() { ss << "2"; });
    op.acquire();
    CHECK(ss.str() == "a = ");
    op.release();
    CHECK(ss.str() == "a = 2");
}

TEST_CASE("Make operations") {
    std::stringstream ss;
    auto greeter = make_operation([&ss]() { ss << "Hello! "; }, [&ss]() { ss << "Goodbye! "; });
    auto hi = make_acquire_operation([&ss]() { ss << "Hi! "; });
    auto bye = make_release_operation([&ss]() { ss << "Bye! "; });

    greeter->acquire();
    hi->acquire();
    bye->acquire();
    CHECK(ss.str() == "Hello! Hi! ");

    greeter->release();
    hi->release();
    bye->release();
    CHECK(ss.str() == "Hello! Hi! Goodbye! Bye! ");
}

struct MyStruct {
    Operation op;
    MyStruct(std::stringstream& ss, int i)
        : op([&ss, i]() { ss << "+" << i; }, [&ss, i]() { ss << "-" << i; }) {}
    double unused{0};
};

TEST_CASE("ForInContainer test") {
    std::stringstream ss;

    // clang-format off
    std::vector<MyStruct> v{{ss, 1}, {ss, 2}, {ss, 3}};
    std::vector<Operation> v2{
        {[&ss]() { ss << "*2"; }, [&ss]() { ss << "/2"; } },
        {[&ss]() { ss << "*3"; }, [&ss]() { ss << "/3"; } }
    };
    ForInContainer<vector<MyStruct>> group(v, [](MyStruct& s) -> Proxy& { return s.op; });
    ForInContainer<vector<Operation>> group2(v2);
    // clang-format on

    group.acquire();
    group2.acquire();
    CHECK(ss.str() == "+1+2+3*2*3");
    group.release();
    group2.release();
    CHECK(ss.str() == "+1+2+3*2*3-1-2-3/2/3");
}

TEST_CASE("make_for_in_container") {
    std::stringstream ss;

    // clang-format off
    std::vector<MyStruct> v{{ss, 1}, {ss, 2}, {ss, 3}};
    std::vector<Operation> v2{
        {[&ss]() { ss << "*2"; }, [&ss]() { ss << "/2"; }},
        {[&ss]() { ss << "*3"; }, [&ss]() { ss << "/3"; }}
    };
    auto group = make_for_in_container(v, [](MyStruct& s) -> Proxy& { return s.op; });
    auto group2 = make_for_in_container(v2);
    // clang-format on

    group->acquire();
    group2->acquire();
    CHECK(ss.str() == "+1+2+3*2*3");
    group->release();
    group2->release();
    CHECK(ss.str() == "+1+2+3*2*3-1-2-3/2/3");
}