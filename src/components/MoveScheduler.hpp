#include <utility>
#include "ChainComponent.hpp"

/* A class that takes a parameter-less lambda and wraps it into
a ChainComponent where the lambda is called on move */
template <class F>
class MoveScheduler : public ChainComponent {
    F f;

  public:
    MoveScheduler(F f) : f(f) {}
    void move(int) final { f(); }
};

template <class F>
auto make_move_scheduler(F&& f) {
    return MoveScheduler<F>(std::forward<F>(f));
}