#pragma once

#include "utils.hpp"

class Proxy {
public:
    virtual void acquire() {}
    virtual void release() {}
};