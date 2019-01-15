#pragma once

#include <type_traits>

template <class... Tags>
struct Context : public Tags... {};

template <class _Context, class _Tag>
struct has_tag : public std::integral_constant<bool, std::is_base_of<_Tag, _Context>::value> {};

template <class _Context, class _Tag>
struct add_tag {
    template <class... Tags>
    static Context<_Tag, Tags...> helper(Context<Tags...>) {}

    using type = decltype(helper(_Context()));
};