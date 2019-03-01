#pragma once

#include <string>
#include "decl_info.hpp"
using std::string;

/*==================================================================================================
  Various utilities to help declare introspection interfaces
==================================================================================================*/

template <class... Tags, class User, class Target, class... Args>
void declare(User& user, string name, Target& target, Args&&... args) {
    user.process_declaration(make_decl_info<Tags...>(target), name, std::forward<Args>(args)...);
}

template <class User, class Provider>
void basic_apply(User& user, Provider& provider) {
    provider.declare_interface(user);
}

namespace decl_utils {
    template <class User, class FilterTag>
    class FilterHelper {
        User& user;

        template <class... Args>
        void filter_dispatch(std::true_type, Args&&... args) {
            user.process_declaration(std::forward<Args>(args)...);
        }

        template <class... Args>
        void filter_dispatch(std::false_type, Args&&... args) {}

      public:
        FilterHelper(User& user) : user(user) {}

        template <class Info, class... Args>
        void process_declaration(Info info, Args&&... args) {
            using has_filter_tag = typename Info::context::template has_tag<FilterTag>;
            filter_dispatch(has_filter_tag(), info, std::forward<Args>(args)...);
        }
    };
}  // namespace decl_utils

template <class Tag, class User, class Provider>
void filter_apply(User& user, Provider& provider) {
    decl_utils::FilterHelper<User, Tag> helper(user);
    basic_apply(helper, provider);
}