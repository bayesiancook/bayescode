#pragma once

#include <type_traits>

namespace is_context {
    struct tag {};
    template <class T>
    struct trait : public std::integral_constant<bool, std::is_base_of<tag, T>::value> {};
}  // namespace is_context

template <class... Tags>
struct Context : is_context::tag, public Tags... {
    using Self = Context<Tags...>;

    template <class Tag>
    struct has_tag : public std::integral_constant<bool, std::is_base_of<Tag, Self>::value> {};

    template <class Tag>
    static Context<Tag, Tags...> add_tag() {
        return Context<Tag, Tags...>();
    }
};


// template <class _Context, class _Context2>
// struct context_union {
//     static_assert(is_context::trait<_Context>::value, "_Context is not a context");

//     static _Context helper(Context<>) {}  // _Context2 is empty

//     template <class Tag, class... Tags>
//     static typename context_union<typename add_tag<_Context, Tag>::type, Context<Tags...>>::type
//     helper(Context<Tag, Tags...>) {}

//     using type = decltype(helper(_Context2()));
// };

// template <class _Context, class _Tag>
// struct remove_tag {
//     static_assert(is_context::trait<_Context>::value, "_Context is not a context");

//     template <class Arg, class... Args>
//     static typename remove_tag<Context<Args...>, _Tag>::type  //
//     select(std::true_type, Context<Arg, Args...>) {}

//     template <class Arg, class... Args>
//     static typename add_tag<typename remove_tag<Context<Args...>, _Tag>::type, Arg>::type  //
//     select(std::false_type, Context<Arg, Args...>) {}

//     static Context<> helper(Context<>) {}  // empty context: nothing to remove

//     template <class Arg, class... Args>
//     static decltype(select(typename std::is_same<Arg, _Tag>::type(), Context<Arg, Args...>()))  //
//     helper(Context<Arg, Args...>) {}

//     using type = decltype(helper(_Context()));
// };