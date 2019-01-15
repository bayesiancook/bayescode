#include "doctest.h"

#include "context.hpp"

struct MyTag {};
struct MyTag2 {};

TEST_CASE("has_tag") {
    Context<MyTag2, MyTag> c;
    Context<MyTag2> d;

    bool c_has_mytag = has_tag<decltype(c), MyTag>::value;
    CHECK(c_has_mytag == true);
    bool d_has_mytag = has_tag<decltype(d), MyTag>::value;
    CHECK(d_has_mytag == false);
}

TEST_CASE("add_tag") {
    Context<> c;

    bool c_has_mytag = has_tag<decltype(c), MyTag>::value;
    CHECK(c_has_mytag == false);

    add_tag<decltype(c), MyTag2>::type d;
    bool d_has_mytag = has_tag<decltype(c), MyTag>::value;
    CHECK(d_has_mytag == false);

    add_tag<decltype(d), MyTag>::type e;
    bool e_has_mytag = has_tag<decltype(e), MyTag>::value;
    CHECK(e_has_mytag == true);
}