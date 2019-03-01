#pragma once

#include <string>
#include "decl_info.hpp"
using std::string;

/*==================================================================================================
  Various utilities to help declare introspection interfaces and apply them
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Basic declare function to be used in declare_interface methods */
template <class... Tags, class User, class Target, class... Args>
void declare(User& user, string name, Target& target, Args&&... args) {
    user.process_declaration(make_decl_info<Tags...>(target), name, std::forward<Args>(args)...);
}

/*--------------------------------------------------------------------------------------------------
  Basic apply: apply interface declaration of provider to user */
template <class User, class Provider>
void basic_apply(User& user, Provider& provider) {
    provider.declare_interface(user);
}

/*--------------------------------------------------------------------------------------------------
  Building bricks to construct application operations */
namespace decl_utils {  // namespace to hide helpers

    template <class User, class Forwarding>
    class Start {
        User& user;

      public:
        Start(User& user) : user(user) {}

        template <class Info, class... Args>
        void process_declaration(Info info, Args&&... args) {
            Forwarding::forward_declaration(user, info, std::forward<Args>(args)...);
        }
    };

    class End {
      public:
        template <class User, class... Args>
        static void forward_declaration(User& user, Args&&... args) {
            user.process_declaration(std::forward<Args>(args)...);
        }
    };

    template <class Tag>
    struct HasTag {
        template <class Info>
        static auto test() {
            static_assert(
                is_decl_info::trait<Info>::value, "Info given to HasTag::test is not a decl info");
            return typename Info::context::template has_tag<Tag>();
        }
    };

    template <class Type>
    struct HasType {
        template <class Info>
        static auto test() {
            static_assert(
                is_decl_info::trait<Info>::value, "Info given to HasType::test is not a decl info");
            return std::is_same<Type, typename Info::target_type>();
        }
    };

    template <class Test, class Forwarding>
    class Filter {
        template <class User, class... Args>
        static void filter_dispatch(std::true_type, User& user, Args&&... args) {
            Forwarding::forward_declaration(user, std::forward<Args>(args)...);
        }

        template <class User, class... Args>
        static void filter_dispatch(std::false_type, User&, Args&&...) {}

      public:
        template <class User, class Info, class... Args>
        static void forward_declaration(User& user, Info info, Args&&... args) {
            static_assert(is_decl_info::trait<Info>::value,
                "Info given to Filter::process_declaration is not a decl info");
            filter_dispatch(Test::template test<Info>(), user, info, std::forward<Args>(args)...);
        }
    };

    template <class Tag, class Forwarding>
    using FilterTag = Filter<HasTag<Tag>, Forwarding>;

    template <class Type, class Forwarding>
    using FilterType = Filter<HasType<Type>, Forwarding>;

}  // namespace decl_utils

/*--------------------------------------------------------------------------------------------------
  Filter apply: allows application of only declarations with a given tag */
template <class Tag, class User, class Provider>
void filter_apply(User& user, Provider& provider) {
    using namespace decl_utils;
    Start<User, FilterTag<Tag, End>> helper(user);
    basic_apply(helper, provider);
}

/*--------------------------------------------------------------------------------------------------
  Typefilter apply: allows application of only declarations with a given target type */
template <class Type, class User, class Provider>
void typefilter_apply(User& user, Provider& provider) {
    decl_utils::Start<User, decl_utils::FilterType<Type, decl_utils::End>> helper(user);
    basic_apply(helper, provider);
}