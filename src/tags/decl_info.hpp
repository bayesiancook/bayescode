#pragma once

#include "context.hpp"
#include "processing_info.hpp"

/*==================================================================================================
  Decl info
  A decl info is data associated to declaration in an introspection interface.
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Decl info trait
  A type trait used to check that a type passed is a Decl info. */
namespace is_decl_info {
    struct tag {};
    template <class T>
    struct trait : public std::integral_constant<bool, std::is_base_of<tag, T>::value> {};
}  // namespace is_decl_info

/*--------------------------------------------------------------------------------------------------
  Decl info class. */
template <class _Context>
struct DeclInfo : is_decl_info::tag {
    using context = _Context;
    static_assert(is_context::trait<context>::value, "context is not a context");

    std::string name;

    DeclInfo(std::string name) : name(name) {}

    template <class Tag>
    auto add_tag() const {  // returns a new object with added tag
        return DeclInfo<decltype(context::template add_tag<Tag>())>(name);
    }

    template <class Tag>
    bool has_tag() const {
        return context::template has_tag<Tag>::value;
    }
};

/*--------------------------------------------------------------------------------------------------
  Decl info factory functions */
template <class... Tags>
auto make_decl_info(std::string name) {  // FIXME: not very useful now that there is no more target
    return DeclInfo<Context<Tags...>>(name);
}
