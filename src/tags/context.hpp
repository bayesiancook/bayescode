#pragma once

#include <type_traits>

template <class Target, class... Tags>
struct Context : public Tags... {
    Target& target;
    Context(Target& target) : target(target) {}
};

template <class _Context, class _Tag>
struct has_tag : public std::integral_constant<bool, std::is_base_of<_Tag, _Context>::value> {};