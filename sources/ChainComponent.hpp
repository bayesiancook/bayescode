#pragma once

class ChainComponent {
public:
    virtual void start() {}
    virtual void move(int) {}
    virtual void after_move(int) {}
    virtual void end() {}
};
