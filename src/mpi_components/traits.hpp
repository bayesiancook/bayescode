#pragma once

#include <vector>

/*
====================================================================================================
  Serializable stuff
==================================================================================================*/
// clang-format off

// is supported by buffers directly (int and double)
template <class T> struct is_default_serializable
    { static bool const value = std::is_same<T, int>::value or std::is_same<T, double>::value; };

// implements template <T> void serialization_interface(T&)
template <class T> struct has_custom_serialization     { static bool const value = false; };

// can be serialized by contiguous chunks (eg, vector<double>)
template <class T> struct is_contiguously_serializable { static bool const value = false; };
// clang-format on

/*
====================================================================================================
  Array stuff
==================================================================================================*/
// clang-format off

// implements size_t size() or T& operator[]()
template <class T> struct has_size            { static bool const value = false; };
template <class T> struct has_access_operator { static bool const value = false; };

// can be partitioned by buffer manager
template <class T> struct is_partitionable 
    { static bool const value = has_size<T>::value and has_access_operator<T>::value; };
// clang-format on

/*
====================================================================================================
  Vector traits
==================================================================================================*/
// clang-format off
template <class T> struct is_contiguously_serializable<std::vector<T>>
    { static bool const value = is_default_serializable<T>::value; };

template <class T> struct has_size<std::vector<T>>            { static bool const value = true; };
template <class T> struct has_access_operator<std::vector<T>> { static bool const value = true; };
// clang-format on