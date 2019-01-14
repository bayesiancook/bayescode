#include "mpi_components/utils.hpp"

using namespace MPI;
using namespace std;

// =========================================================================
//   GLOBAL DEFS
// =========================================================================

// forwarding based on target traits
namespace needs_name {
    // clang-format off
    struct tag {};

    template <class T> struct trait :
        public integral_constant<bool, is_base_of<tag, T>::value> {};
    // clang-format on

    template <class Target, class... Args>
    void continuation(Target& target, Args&&... args) {
        target.add(std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void select(true_type, Target& target, string name, Args&&... args) {
        continuation(target, name, std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void select(false_type, Target& target, string, Args&&... args) {
        continuation(target, std::forward<Args>(args)...);
    }
}  // namespace needs_name

template <class Target, class... Args>
void declare(Target& target, string name, Args&&... args) {
    needs_name::select(needs_name::trait<Target>(), target, name, std::forward<Args>(args)...);
}


// =========================================================================
//   TOOLS
// =========================================================================

struct MyPrinter /*: needs_name::tag*/ {
    void add(string name, double value) {
        p->message("Variable %s has value %f", name.c_str(), value);
    }

    void add(string name, int value) {
        p->message("Variable %s has value %d", name.c_str(), value);
    }
};

// clang-format off
template <> struct needs_name::trait<MyPrinter> : public true_type {};
// clang-format on


// =========================================================================
//   MODEL
// =========================================================================

class MyModel {
    double x = 2.3;
    int i = 5;

  public:
    template <class Target>
    void interface(Target& target) {
        declare(target, "x", x);
        declare(target, "i", i);
    }
};


void compute(int, char**) {
    MyModel m;
    MyPrinter p;
    m.interface(p);
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }