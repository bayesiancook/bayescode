#pragma once

/*==================================================================================================
  Building bricks to construct application operations
==================================================================================================*/
namespace processing {  // namespace to hide helpers

    /*----------------------------------------------------------------------------------------------
      Ending bricks. Objects that call the final user. */
    class End {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(
                declinfo.name, declinfo.target, std::forward<Args>(args)...);
        }
    };

    class FullNameEnd {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(
                prinfo.name + declinfo.name, declinfo.target, std::forward<Args>(args)...);
        }
    };

    class NoNameEnd {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(declinfo.target, std::forward<Args>(args)...);
        }
    };

    /*----------------------------------------------------------------------------------------------
      Tests: boolean properties on types. */
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

    /*----------------------------------------------------------------------------------------------
      Processing bricks. */
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

    template <class Test, class Forwarding, bool recursive>
    class Unroll {
        template <class PrInfo, class DeclInfo, class... Args>  // to be unrolled once
        static void filter_dispatch(
            std::true_type, std::false_type, PrInfo prinfo, DeclInfo declinfo, Args&&...) {
            // TODO: fix redundant info regarding current processing (prinfo + current class) ?
            prinfo.name += declinfo.name + "_";
            declinfo.target.declare_interface(make_processing_info<Forwarding>(prinfo.user));
        }

        template <class PrInfo, class DeclInfo, class... Args>  // to be unrolled
        static void filter_dispatch(
            std::true_type, std::true_type, PrInfo prinfo, DeclInfo declinfo, Args&&...) {
            // TODO: fix redundant info regarding current processing (prinfo + current class) ?
            prinfo.name += declinfo.name + "_";
            declinfo.target.declare_interface(prinfo);
        }

        template <class Anything, class... Args>  // not to be unrolled
        static void filter_dispatch(std::false_type, Anything, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO check types
            filter_dispatch(Test::template test<DeclInfo>(),
                std::integral_constant<bool, recursive>(), prinfo, declinfo,
                std::forward<Args>(args)...);
        }
    };

    template <class Test, class Forwarding>
    using SimpleUnroll = Unroll<Test, Forwarding, false>;

    template <class Test, class Forwarding>
    using RecursiveUnroll = Unroll<Test, Forwarding, true>;

    template <class Tag, class Forwarding>
    using UnrollIf = SimpleUnroll<HasTag<Tag>, Forwarding>;

}  // namespace processing
