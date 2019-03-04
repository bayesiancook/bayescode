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

    class NoNameEnd {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(declinfo.target, std::forward<Args>(args)...);
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

    template <class Test, class Forwarding>
    class SimpleUnroll {
        template <class PrInfo, class DeclInfo, class... Args>  // to be unrolled
        static void filter_dispatch(std::true_type, PrInfo prinfo, DeclInfo declinfo, Args&&...) {
            // NOTE: name and args are discarded! (FIXME?)
            // TODO: fix redundant info regarding current processing (prinfo + current class)
            declinfo.target.declare_interface(make_processing_info<Forwarding>(prinfo.user));
        }

        template <class... Args>  // not to be unrolled
        static void filter_dispatch(std::false_type, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO check types
            filter_dispatch(
                Test::template test<DeclInfo>(), prinfo, declinfo, std::forward<Args>(args)...);
        }
    };

    template <class Test, class Forwarding>
    class RecursiveUnroll {
        template <class PrInfo, class DeclInfo, class... Args>  // to be unrolled
        static void filter_dispatch(std::true_type, PrInfo prinfo, DeclInfo declinfo, Args&&...) {
            // NOTE: name and args are discarded! (FIXME?)
            // TODO: fix redundant info regarding current processing (prinfo + current class)
            declinfo.target.declare_interface(prinfo);
        }

        template <class... Args>  // not to be unrolled
        static void filter_dispatch(std::false_type, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO check types
            filter_dispatch(
                Test::template test<DeclInfo>(), prinfo, declinfo, std::forward<Args>(args)...);
        }
    };


    template <class Tag, class Forwarding>
    using UnrollIf = SimpleUnroll<HasTag<Tag>, Forwarding>;

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

/*--------------------------------------------------------------------------------------------------
  Single unroll of structures with tag Tag */
template <class Tag, class User, class Provider>
void unrollif_apply(User& user, Provider& provider) {
    using namespace decl_utils;
    auto processing_info = make_processing_info<UnrollIf<Tag, End>>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Recursive unroll of structures with tag Tag */
template <class Tag, class User, class Provider>
void recif_apply(User& user, Provider& provider) {
    using namespace decl_utils;
    auto processing_info = make_processing_info<RecursiveUnroll<HasTag<Tag>, End>>(user);
    provider.declare_interface(processing_info);
}