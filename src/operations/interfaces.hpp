#pragma once

class Proxy {
  public:
    virtual void acquire() {}
    virtual void release() {}
    virtual ~Proxy() = default;
};