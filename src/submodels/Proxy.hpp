#pragma once

template <class T>
class Proxy {
    virtual T _get() = 0;

  protected:
    // protected non-virtual destructor, as this interface is not
    // meant to be sued with owning pointers
    ~Proxy() = default;

  public:
    virtual void gather() = 0;

    T get() {
#ifndef NDEBUG
        auto tmp = _get();
        gather();
        assert(_get() == tmp);
#endif
        return _get();
    }
};