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
            prinfo.user.process_declaration(declinfo.name, std::forward<Args>(args)...);
        }
    };

    class FullNameEnd {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(
                prinfo.name + declinfo.name, std::forward<Args>(args)...);
        }
    };

    class NoNameEnd {
      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            prinfo.user.process_declaration(std::forward<Args>(args)...);
        }
    };

    /*----------------------------------------------------------------------------------------------
      Trait mappers. */
    template <class Tag>
    struct HasTag {
        template <class DeclInfo, class... Args>
        using constant = typename DeclInfo::context::template has_tag<Tag>;
    };

    template <class Type>
    struct HasType {
        template <class DeclInfo, class Target, class... Args>
        using constant = std::is_same<Type, std::remove_reference_t<Target>>;
    };

    template <template <class> class Trait>
    struct HasTrait {
        template <class DeclInfo, class Target, class... Args>
        using constant = Trait<std::remove_reference_t<Target>>;
    };

    struct True {
        template <class... Args>
        using constant = std::true_type;
    };

    template <class Mapper>
    struct Not {
        template <class DeclInfo, class Target, class... Args>
        using constant = std::integral_constant<bool,
            not Mapper::template constant<DeclInfo, Target, Args...>::value>;
    };

    template <class TraitMapper, class... Args>
    auto get_constant_value() {
        return typename TraitMapper::template constant<Args...>();
    }

    /*----------------------------------------------------------------------------------------------
      Processing bricks. */
    template <class TraitMapper, class Forwarding>
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
            filter_dispatch(get_constant_value<TraitMapper, DeclInfo, Args...>(), prinfo, declinfo,
                std::forward<Args>(args)...);
        }
    };

    template <class Tag, class Forwarding>
    using FilterTag = Filter<HasTag<Tag>, Forwarding>;

    template <class Type, class Forwarding>
    using FilterType = Filter<HasType<Type>, Forwarding>;

    template <class TraitMapper, class Forwarding, bool recursive>
    class Unroll {
        template <class PrInfo, class DeclInfo, class Target, class... Args>  // to be unrolled once
        static void filter_dispatch(std::true_type, std::false_type, PrInfo prinfo,
            DeclInfo declinfo, Target& target, Args&&...) {
            // TODO: fix redundant info regarding current processing (prinfo + current class) ?
            prinfo.name += declinfo.name + "_";
            target.declare_interface(make_processing_info<Forwarding>(prinfo.user));
        }

        template <class PrInfo, class DeclInfo, class Target, class... Args>  // to be unrolled
        static void filter_dispatch(std::true_type, std::true_type, PrInfo prinfo,
            DeclInfo declinfo, Target& target, Args&&...) {
            // TODO: fix redundant info regarding current processing (prinfo + current class) ?
            auto new_prinfo = make_processing_info<Unroll<TraitMapper, Forwarding, recursive>>(
                prinfo.user, prinfo.name + declinfo.name + "_");
            target.declare_interface(new_prinfo);
        }

        template <class Anything, class... Args>  // not to be unrolled
        static void filter_dispatch(std::false_type, Anything, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO check types
            filter_dispatch(get_constant_value<TraitMapper, DeclInfo, Args...>(),
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
