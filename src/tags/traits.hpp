#pragma once

#include <type_traits>
#include "global/logging.hpp"

template <class, class>
struct ProcessingInfo;
struct Ignore;

/*--------------------------------------------------------------------------------------------------
  Trait to check that a type has an interface */
template <class T>
class has_interface_helper {
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