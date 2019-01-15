#include "doctest.h"

#include "context.hpp"

TEST_CASE("has_tag") {
    struct MyTag {};
    struct MyTag2 {};
    double a = 2;
    Context<double, MyTag2, MyTag> c{a};
    Context<double, MyTag2> d{a};

    bool c_has_mytag = has_tag<decltype(c), MyTag>::value;
    CHECK(c_has_mytag == true);
    bool d_has_mytag = has_tag<decltype(d), MyTag>::value;
    CHECK(d_has_mytag == false);
}