#pragma once

#include <type_traits>

/*==================================================================================================
  Context
  A context is a set of type tags.
  It is meant to be used to encode metadata for individual declarations in introspection interfaces.
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Context trait.
  A type trait used to check that a type passed is a Context. */
namespace is_context {
    struct tag {};
    template <class T>
    struct trait : public std::integral_constant<bool, std::is_base_of<tag, T>::value> {};
}  // namespace is_context

/*--------------------------------------------------------------------------------------------------
  Context class. */
template <class... Tags>
class Context : is_context::tag, public Tags... {
    using Self = Context<Tags...>;

    template <class... OtherTags>
    static auto merge_helper(Context<OtherTags...>) {
        return Context<Tags..., OtherTags...>();
    }

  public:
    template <class Tag>
    struct has_tag : public std::integral_constant<bool, std::is_base_of<Tag, Self>::value> {};

    template <class Tag>
    static auto add_tag() {
        return Context<Tag, Tags...>();
    }

    template <class OtherContext>
    static auto merge(OtherContext) {
        static_assert(is_context::trait<OtherContext>::value, "OtherContext is not a context");
        return merge_helper(OtherContext());
    }
};
