#pragma once

#include <type_traits>
#include "global/logging.hpp"

template <class, class>
struct ProcessingInfo;
struct Ignore;

/*--------------------------------------------------------------------------------------------------
  Trait to check that a type has an interface */
template <class T>
class has_interface_helper {  // @todo: rewrite this trait properly!
    template <class T2>
    static constexpr
        typename std::is_same<void, decltype(std::declval<T2>().declare_interface(
                                        std::declval<ProcessingInfo<int, Ignore>>()))>::type
        helper(int) {
        return std::true_type();
    }

    template <class>
    static constexpr auto helper(float) {
        return std::false_type();
    }

  public:
    using constant = decltype(helper<T>(0));
};

template <class T>
using has_interface = typename has_interface_helper<T>::constant;

template <class T>
struct external_interface;

template <class... T>
using my_void_t = void;

template <class T, class = void>
struct has_external_interface : std::false_type {};

template <class T>
struct has_external_interface<T, my_void_t<decltype(external_interface<T>{})>> : std::true_type {};

template <class T>
struct has_either_interface
    : std::integral_constant<bool, has_interface<T>::value || has_external_interface<T>::value> {};