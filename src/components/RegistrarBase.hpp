#pragma once

#include <assert.h>
#include <set>
#include <string>

class Partition;

using filter_t = std::set<std::string>;

// TODO move into a logger class somewhere
#ifdef __GNUG__
#include <cxxabi.h>
static inline std::string demangle(const char* name) {
    int status{0};
    std::unique_ptr<char, void (*)(void*)> res{
        abi::__cxa_demangle(name, NULL, NULL, &status), std::free};
    return (status == 0) ? res.get() : name;
}
#else
static std::string demangle(const char* name) { return name; }
#endif

template <class T>
class RegistrarBase {
    mutable filter_t _filter;

  public:
    template <class... Args>
    void register_element(std::string s, Args... args) {
        std::stringstream ss;
        ss << "\e[1m\e[31mError\e[0m| \e[33mregister_element\e[0m overload for element \e[32m" << s
           << "\e[0m not implemented for class \e[32m" << demangle(typeid(T).name())
           << "\e[0m\n     \\ overload parameters are: \e[32m";
        std::vector<std::string> types{demangle(typeid(Args).name())...};
        for (auto t : types) { ss << t << " "; }
        ss << "\e[0m\n";
        std::cerr << ss.str();
        exit(1);
    }

    RegistrarBase() = default;

    template <class... Args>
    void add(std::string name, Args&&... args) {
        if ((_filter.size() == 0) or (_filter.find(name) != _filter.end())) {
            static_cast<T*>(this)->register_element(name, std::forward<Args>(args)...);
        }
    }

    template <class M>
    void register_from_method(
        M& ref, void (M::*f)(RegistrarBase<T>&), filter_t filter = filter_t{}) {
        /* -- */
        assert(filter.size() > 0);
        _filter = filter;  // local copy
        (ref.*f)(*this);
        _filter.clear();
    }
};
