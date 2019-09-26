#pragma once

template <class T, class... Args>
class Proxy {
    virtual T _get(Args... args) = 0;

  protected:
    // protected non-virtual destructor, as this interface is not
    // meant to be used with owning pointers
    ~Proxy() = default;

  public:
    // @fixme? maybe this should be a private virtual _gather for consistency
    virtual void gather() = 0;

    T get(Args... args) {
#ifndef NDEBUG
        auto tmp = _get(args...);
        gather();
        assert(_get(args...) == tmp);
#endif
        return _get(args...);
    }
};