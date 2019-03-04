#pragma once

#include <string>
#include "decl_info.hpp"
using std::string;

/*==================================================================================================
  Various utilities to help declare introspection interfaces and apply them
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Basic declare function to be used in declare_interface methods */
template <class... Tags, class PrInfo, class Target, class... Args>
void declare(PrInfo& processing_info, string name, Target& target, Args&&... args) {
    auto decl_info = make_decl_info<Tags...>(target, name);
    PrInfo::processing::forward_declaration(
        processing_info, decl_info, std::forward<Args>(args)...);
}

/*--------------------------------------------------------------------------------------------------
  Building bricks to construct application operations */
namespace decl_utils {  // namespace to hide helpers

    class End {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(
                declinfo.name, declinfo.target, std::forward<Args>(args)...);
        }
    };

    //     class NoNameEnd {
    //       public:
    //         template <class User, class Info, class... Args>
    //         static void forward_declaration(User& user, Info info, string, Args&&... args) {
    //             user.process_declaration(info.target, std::forward<Args>(args)...);
    //         }
    //     };

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
        template <class... Args>
        static void filter_dispatch(std::true_type, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

        template <class... Args>
        static void filter_dispatch(std::false_type, Args&&...) {}

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO: check infos are infos
            filter_dispatch(
                Test::template test<DeclInfo>(), prinfo, declinfo, std::forward<Args>(args)...);
        }
    };

    template <class Tag, class Forwarding>
    using FilterTag = Filter<HasTag<Tag>, Forwarding>;

    template <class Type, class Forwarding>
    using FilterType = Filter<HasType<Type>, Forwarding>;

    //     template <class Test, class Forwarding>
    //     class SimpleUnroll {
    //         template <class User, class Info, class... Args>  // to be unrolled
    //         static void filter_dispatch(std::true_type, User& user, Info info, Args&&...) {
    //             // NOTE: name and args are discarded! (FIXME?)
    //             info.target.template declare_interface<Forwarding>(user);
    //         }

    //         template <class User, class... Args>  // not to be unrolled
    //         static void filter_dispatch(std::false_type, User& user, Args&&... args) {
    //             Forwarding::forward_declaration(user, std::forward<Args>(args)...);
    //         }

    //       public:
    //         template <class User, class Info, class... Args>
    //         static void forward_declaration(User& user, Info info, Args&&... args) {
    //             static_assert(is_decl_info::trait<Info>::value,
    //                 "Info given to Unroll::process_declaration is not a decl info");
    //             filter_dispatch(Test::template test<Info>(), user, info,
    //             std::forward<Args>(args)...);
    //         }
    //     };

    //     template <class Test, class Forwarding>
    //     class RecursiveUnroll {
    //         template <class User, class Info, class... Args>  // to be unrolled
    //         static void filter_dispatch(std::true_type, User& user, Info info, Args&&...) {
    //             // NOTE: name and args are discarded! (FIXME?)
    //             info.target.template declare_interface<RecursiveUnroll<Test, Forwarding>>(user);
    //         }

    //         template <class User, class... Args>  // not to be unrolled
    //         static void filter_dispatch(std::false_type, User& user, Args&&... args) {
    //             Forwarding::forward_declaration(user, std::forward<Args>(args)...);
    //         }

    //       public:
    //         template <class User, class Info, class... Args>
    //         static void forward_declaration(User& user, Info info, Args&&... args) {
    //             static_assert(is_decl_info::trait<Info>::value,
    //                 "Info given to Unroll::process_declaration is not a decl info");
    //             filter_dispatch(Test::template test<Info>(), user, info,
    //             std::forward<Args>(args)...);
    //         }
    //     };

    //     template <class Tag, class Forwarding>
    //     using UnrollIf = SimpleUnroll<HasTag<Tag>, Forwarding>;

}  // namespace decl_utils

/*--------------------------------------------------------------------------------------------------
  Basic apply: apply interface declaration of provider to user */
template <class User, class Provider>
void basic_apply(User& user, Provider& provider) {
    auto processing_info = make_processing_info<decl_utils::End>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Filter apply: allows application of only declarations with a given tag */
template <class Tag, class User, class Provider>
void filter_apply(User& user, Provider& provider) {
    using namespace decl_utils;
    auto processing_info = make_processing_info<FilterTag<Tag, End>>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Typefilter apply: allows application of only declarations with a given target type */
template <class Type, class User, class Provider>
void typefilter_apply(User& user, Provider& provider) {
    using namespace decl_utils;
    auto processing_info = make_processing_info<FilterType<Type, End>>(user);
    provider.declare_interface(processing_info);
}

// /*--------------------------------------------------------------------------------------------------
//   Single unroll of structures with tag Tag */
// template <class Tag, class User, class Provider>
// void unrollif_apply(User& user, Provider& provider) {
//     using namespace decl_utils;
//     provider.template declare_interface<UnrollIf<Tag, End>>(user);
// }

// /*--------------------------------------------------------------------------------------------------
//   Recursive unroll of structures with tag Tag */
// template <class Tag, class User, class Provider>
// void recif_apply(User& user, Provider& provider) {
//     using namespace decl_utils;
//     provider.template declare_interface<RecursiveUnroll<HasTag<Tag>, End>>(user);
// }