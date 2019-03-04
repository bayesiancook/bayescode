#include "doctest.h"

#include "decl_utils.hpp"

using namespace std;

struct MyTag {};
struct MyTag2 {};
struct MyTag3 {};
struct MyTag4 {};

/*==================================================================================================
  Contexts
==================================================================================================*/
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

/*==================================================================================================
  DeclInfo
==================================================================================================*/
TEST_CASE("DeclInfo basic usage") {
    double a = 2.33;
    auto i = make_decl_info<MyTag, MyTag2>(a, "a");

    CHECK(i.target == 2.33);
    bool target_is_double = is_same<double, decltype(i)::target_type>::value;
    CHECK(target_is_double == true);
    bool context_is_tag1_tag2 = is_same<Context<MyTag, MyTag2>, decltype(i)::context>::value;
    CHECK(context_is_tag1_tag2 == true);
}

TEST_CASE("DeclInfo::add_tag") {
    double a = 2.33;
    auto i = make_decl_info<MyTag>(a, "a");

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

TEST_CASE("DeclInfo::has_tag") {
    double a = 3.22;
    auto i = make_decl_info<MyTag2>(a, "a");

    CHECK(!i.has_tag<MyTag>());
    CHECK(i.has_tag<MyTag2>());
    CHECK(!i.has_tag<MyTag3>());
}

/*==================================================================================================
  Declaration and use
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Basic cases */
struct Provider {
    int a{1}, b{2}, c{3};

    template <class Info>
    void declare_interface(Info info) {
        declare<MyTag2>(info, "a", a);
        declare<MyTag>(info, "b", b);
        declare<MyTag, MyTag2>(info, "c", c);
    }
};

struct User {
    int sum{0};

    void process_declaration(string, int value) { sum += value; }
};

TEST_CASE("Decl utils base") {
    Provider p;
    User u;

    basic_apply(u, p);
    CHECK(u.sum == 6);
}

/*--------------------------------------------------------------------------------------------------
  Filtering */
TEST_CASE("Decl utils filter") {
    Provider p;
    User u;

    filter_apply<MyTag>(u, p);
    CHECK(u.sum == 5);
    u.sum = 0;

    filter_apply<MyTag2>(u, p);
    CHECK(u.sum == 4);
    u.sum = 0;

    filter_apply<MyTag3>(u, p);
    CHECK(u.sum == 0);
}

struct Provider2 {
    int a{2}, b{13};
    string c;

    template <class Info>
    void declare_interface(Info info) {
        declare<MyTag>(info, "a", a);
        declare<MyTag>(info, "b", b);
        declare<MyTag2>(info, "c", c);
    }
};

TEST_CASE("Filter apply: check that options that would not compile are not compiled if filtered") {
    Provider2 p;
    User u;

    filter_apply<MyTag>(u, p);
    CHECK(u.sum == 15);
}

TEST_CASE("Filter by type") {
    Provider2 p;
    User u;
    typefilter_apply<int>(u, p);

    CHECK(u.sum == 15);
}

// /*--------------------------------------------------------------------------------------------------
//   Forwarding of other arguments */
// struct Provider3 {
//     int a{2}, b{5}, c{11};

//     template <class Info>
//     void declare_interface(Info info) {
//         declare(info, "a", a, true);
//         declare(info, "b", b, false);
//         declare(info, "c", c, true);
//     }
// };

// struct User2 {
//     int sum{0};

//     void process_declaration(string, int value, bool toggle) {
//         if (toggle) { sum += value; }
//     }
// };

// TEST_CASE("Argument forwarding") {
//     User2 u;
//     Provider3 p;

//     basic_apply(u, p);
//     CHECK(u.sum == 13);
// }

// /*--------------------------------------------------------------------------------------------------
//   Recursive structures */
// struct Recursive {};

// struct ProviderRec {
//     Provider p;
//     int a{213};

//     template <class Info>
//     void declare_interface(Info info) {
//         declare<Recursive>(info, "p", p);
//         declare(info, "a", a);
//     }
// };

// TEST_CASE("Unroll one level") {
//     User u;
//     ProviderRec p;

//     unrollif_apply<Recursive>(u, p);
//     CHECK(u.sum == 219);
// }

// struct ProviderRec2 {
//     ProviderRec p;
//     int a{15}, b{3};

//     template <class Info>
//     void declare_interface(Info info) {
//         declare<Recursive>(info, "p", p);
//         declare(info, "a", a);
//         declare(info, "b", b);
//     }
// };

// TEST_CASE("Recursive unroll") {
//     User u, u2;
//     ProviderRec2 p;

//     // single unroll
//     using namespace decl_utils;
//     p.declare_interface<SimpleUnroll<HasTag<Recursive>, Filter<HasType<int>, End>>>(u);
//     CHECK(u.sum == 231);

//     // recursive unroll
//     recif_apply<Recursive>(u2, p);
//     CHECK(u2.sum == 237);
// }

/*--------------------------------------------------------------------------------------------------
  No name */
struct User3 {
    int sum{0};

    void process_declaration(int value) { sum += value; }
};

TEST_CASE("NoName end brick") {
    User3 u;
    Provider p;

    p.declare_interface(make_processing_info<decl_utils::NoNameEnd>(u));
}