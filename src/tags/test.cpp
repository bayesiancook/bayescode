#include "doctest.h"

#include "decl_info.hpp"
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

    bool c_has_mytag = decltype(c)::has_tag<MyTag>::value;
    CHECK(c_has_mytag == true);
    bool d_has_mytag = decltype(d)::has_tag<MyTag>::value;
    CHECK(d_has_mytag == false);
}

TEST_CASE("add_tag") {
    Context<> c;

    bool c_has_mytag = decltype(c)::has_tag<MyTag>::value;
    CHECK(c_has_mytag == false);

    auto d = c.add_tag<MyTag2>();
    bool d_has_mytag = decltype(c)::has_tag<MyTag>::value;
    CHECK(d_has_mytag == false);

    auto e = d.add_tag<MyTag>();
    bool e_has_mytag = decltype(e)::has_tag<MyTag>::value;
    CHECK(e_has_mytag == true);
}

TEST_CASE("context merge") {
    Context<MyTag> c;
    Context<MyTag2, MyTag3> d;
    auto e = c.merge(d);

    bool e_has_mytag = decltype(e)::has_tag<MyTag>::value;
    bool e_has_mytag2 = decltype(e)::has_tag<MyTag2>::value;
    bool e_has_mytag3 = decltype(e)::has_tag<MyTag3>::value;
    bool e_has_mytag4 = decltype(e)::has_tag<MyTag4>::value;
    CHECK(e_has_mytag == true);
    CHECK(e_has_mytag2 == true);
    CHECK(e_has_mytag3 == true);
    CHECK(e_has_mytag4 == false);
}

// TEST_CASE("remove_tag") {
//     Context<MyTag> c;
//     bool c_has_mytag = decltype(c)::has_tag<MyTag>::value;
//     bool c_has_mytag2 = decltype(c)::has_tag<MyTag2>::value;
//     bool c_has_mytag3 = decltype(c)::has_tag<MyTag3>::value;
//     CHECK(c_has_mytag == true);
//     CHECK(c_has_mytag2 == true);
//     CHECK(c_has_mytag3 == false);

//     auto d = c.remove_tag<MyTag>();
//     bool d_has_mytag = decltype(d)::has_tag<MyTag>::value;
//     bool d_has_mytag2 = decltype(d)::has_tag<MyTag2>::value;
//     bool d_has_mytag3 = decltype(d)::has_tag<MyTag3>::value;
//     CHECK(d_has_mytag == false);
//     CHECK(d_has_mytag2 == true);
//     CHECK(d_has_mytag3 == false);

//     auto e = d.remove_tag<MyTag3>();
//     bool e_has_mytag = decltype(e)::has_tag<MyTag>::value;
//     bool e_has_mytag2 = decltype(e)::has_tag<MyTag2>::value;
//     bool e_has_mytag3 = decltype(e)::has_tag<MyTag3>::value;
//     CHECK(e_has_mytag == false);
//     CHECK(e_has_mytag2 == true);
//     CHECK(e_has_mytag3 == false);
// }

TEST_CASE("DeclInfo basic usage") {
    double a = 2.33;
    auto i = make_decl_info<MyTag, MyTag2>(a);

    CHECK(i.target == 2.33);
    bool target_is_double = is_same<double, decltype(i)::target_type>::value;
    CHECK(target_is_double == true);
    bool context_is_tag1_tag2 = is_same<Context<MyTag, MyTag2>, decltype(i)::context>::value;
    CHECK(context_is_tag1_tag2 == true);
}

TEST_CASE("DeclInfo::add_tag") {
    double a = 2.33;
    auto i = make_decl_info<MyTag>(a);

    CHECK(i.target == 2.33);
    bool i_has_mytag = decltype(i)::context::has_tag<MyTag>::value;
    bool i_has_mytag2 = decltype(i)::context::has_tag<MyTag2>::value;
    CHECK(i_has_mytag == true);
    CHECK(i_has_mytag2 == false);

    auto j = i.add_tag<MyTag2>();
    CHECK(j.target == 2.33);
    bool j_has_mytag = decltype(j)::context::has_tag<MyTag>::value;
    bool j_has_mytag2 = decltype(j)::context::has_tag<MyTag2>::value;
    CHECK(j_has_mytag == true);
    CHECK(j_has_mytag2 == true);
}