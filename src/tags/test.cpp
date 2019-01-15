#include "doctest.h"

#include "context.hpp"
// #include <typeinfo>
// #include <iostream>

using namespace std;

struct MyTag {};
struct MyTag2 {};
struct MyTag3 {};
struct MyTag4 {};

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

TEST_CASE("context_union") {
    Context<MyTag> c;
    Context<MyTag2, MyTag3> d;
    context_union<decltype(c), decltype(d)>::type e;
    // cout << typeid(e).name() << endl;

    bool e_has_mytag = has_tag<decltype(e), MyTag>::value;
    bool e_has_mytag2 = has_tag<decltype(e), MyTag2>::value;
    bool e_has_mytag3 = has_tag<decltype(e), MyTag3>::value;
    bool e_has_mytag4 = has_tag<decltype(e), MyTag4>::value;
    CHECK(e_has_mytag == true);
    CHECK(e_has_mytag2 == true);
    CHECK(e_has_mytag3 == true);
    CHECK(e_has_mytag4 == false);
}