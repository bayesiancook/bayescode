#pragma once

#include "traits.hpp"

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
            // INFO("Forwarding declaration {}", prinfo.name + declinfo.name);
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

    template <class Mapper1, class Mapper2>
    struct Or {
        template <class DeclInfo, class Target, class... Args>
        using constant = std::integral_constant<bool,
            Mapper1::template constant<DeclInfo, Target, Args...>::value or
                Mapper2::template constant<DeclInfo, Target, Args...>::value>;
    };

    template <class Mapper1, class Mapper2>
    struct And {
        template <class DeclInfo, class Target, class... Args>
        using constant = std::integral_constant<bool,
            Mapper1::template constant<DeclInfo, Target, Args...>::value and
                Mapper2::template constant<DeclInfo, Target, Args...>::value>;
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

    /*----------------------------------------------------------------------------------------------
      Function to select and run interface declaration */
    // TODO: move to other file?
    template <class Prinfo, class Target>  // has_interface
    void call_interface_helper(std::true_type, Prinfo prinfo, Target& target) {
        target.declare_interface(prinfo);
    }

    template <class Prinfo, class Target>  // does not have interface (then call external)
    void call_interface_helper(std::false_type, Prinfo prinfo, Target& target) {
        external_interface<Target>::declare_interface(prinfo, target);
    }

    template <class Prinfo, class Target>
    void call_interface(Prinfo prinfo, Target& target) {
        call_interface_helper(has_interface<Target>(), prinfo, target);
    }

    template <class TraitMapper, class Forwarding, bool recursive>
    class Unroll {
        template <class PrInfo, class DeclInfo, class Target, class... Args>  // to be unrolled once
        static void dispatch(std::true_type, std::false_type, PrInfo prinfo, DeclInfo declinfo,
            Target& target, Args&&...) {
            // TODO: fix redundant info regarding current processing (prinfo + current class) ?
            auto new_prinfo = make_processing_info<Forwarding>(prinfo.user, declinfo.name + "_");
            call_interface(new_prinfo, target);
        }

        template <class PrInfo, class DeclInfo, class Target, class... Args>  // to be unrolled
        static void dispatch(std::true_type, std::true_type, PrInfo prinfo, DeclInfo declinfo,
            Target& target, Args&&...) {
            /* -- */
            DEBUG("Declaration {} must be unrolled. Context is {}", declinfo.name, prinfo.name);
            auto new_prinfo = make_processing_info<Unroll<TraitMapper, Forwarding, recursive>>(
                prinfo.user, prinfo.name + declinfo.name + "_");
            call_interface(new_prinfo, target);
        }

        template <class Anything, class... Args>  // not to be unrolled
        static void dispatch(std::false_type, Anything, Args&&... args) {
            Forwarding::forward_declaration(std::forward<Args>(args)...);
        }

      public:
        template <class PrInfo, class DeclInfo, class... Args>
        static void forward_declaration(PrInfo prinfo, DeclInfo declinfo, Args&&... args) {
            // TODO check types
            dispatch(get_constant_value<TraitMapper, DeclInfo, Args...>(),
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
