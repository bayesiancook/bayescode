#include "doctest.h"

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
    group.add(std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new MyProxy(b))));

    group.acquire();  // 3 b and 2 a
    group.release();
    CHECK(a == 24);
    CHECK(b == 120);
}