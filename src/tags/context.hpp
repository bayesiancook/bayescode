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

template <class _Context, class _Context2>
struct context_union {
    static _Context helper(Context<>) {}  // _Context2 is empty

    template <class Tag, class... Tags>
    static typename context_union<typename add_tag<_Context, Tag>::type, Context<Tags...>>::type
    helper(Context<Tag, Tags...>) {}

    using type = decltype(helper(_Context2()));
};

// template <class _Context, class _Tag>
// struct remove_tag {
//     static _Context helper() {}  // empty context: nothing to remove

//     template <class Arg, class... Args>
//     static decltype(select(std::is_same<Arg, _Tag>(), arg, args)) helper(Arg arg, Args... args)
//     {}

//     template <class Arg, class... Args>
//     static  select(std::true_type, Arg, Args...) {}
// };