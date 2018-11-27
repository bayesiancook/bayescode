#include "doctest.h"

#include "proxies.hpp"

using namespace std;

int a = 1;
int b = 1;

class MyProxy : public Proxy {
    int& target;

  public:
    MyProxy(int& target) : target(target) {}
    void acquire() final { target += 1; }
    void release() final { target *= 2; }
};

TEST_CASE("ForAll test") {
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