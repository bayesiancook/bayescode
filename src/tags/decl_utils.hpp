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
    template <class Tag>
    struct HasTag {
        template <class Info>
        static auto test() {
            return typename Info::context::template has_tag<Tag>();
        }
    };

    template <class Type>
    struct HasType {
        template <class Info>
        static auto test() {
            return std::is_same<Type, typename Info::target_type>();
        }
    };

    template <class User, class Test>
    class Filter {
        User& user;

        template <class... Args>
        void filter_dispatch(std::true_type, Args&&... args) {
            user.process_declaration(std::forward<Args>(args)...);
        }

        template <class... Args>
        void filter_dispatch(std::false_type, Args&&... args) {}

      public:
        Filter(User& user) : user(user) {}

        template <class Info, class... Args>
        void process_declaration(Info info, Args&&... args) {
            filter_dispatch(Test::template test<Info>(), info, std::forward<Args>(args)...);
        }
    };

    template <class User, class Tag>
    using FilterTag = Filter<User, HasTag<Tag>>;

    template <class User, class Type>
    using FilterType = Filter<User, HasType<Type>>;
}  // namespace decl_utils

/*--------------------------------------------------------------------------------------------------
  Filter apply: allows application of only declarations with a given tag */
template <class Tag, class User, class Provider>
void filter_apply(User& user, Provider& provider) {
    decl_utils::FilterTag<User, Tag> helper(user);
    basic_apply(helper, provider);
}

/*--------------------------------------------------------------------------------------------------
  Typefilter apply: allows application of only declarations with a given target type */
template <class Type, class User, class Provider>
void typefilter_apply(User& user, Provider& provider) {
    decl_utils::FilterType<User, Type> helper(user);
    basic_apply(helper, provider);
}