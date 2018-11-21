#pragma once

class ChainComponent {
  public:
    virtual void start() {}
    virtual void move(int) {}
    virtual void savepoint(int) {}
    virtual void end() {}
    virtual ~ChainComponent() = default;
};
