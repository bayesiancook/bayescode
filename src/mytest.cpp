#include "mpi_components/utils.hpp"

using namespace MPI;
using namespace std;

// =========================================================================
//   GLOBAL DEFS
// =========================================================================

template <class Target, class... Args>
void add_to_target(Target& target, Args&&... args) {
    target.add(std::forward<Args>(args)...);
}

namespace needs_name {
    // clang-format off
    struct tag {};
    template <class T> struct trait :
        public integral_constant<bool, is_base_of<tag, T>::value> {};
    // clang-format on

    template <class... Args>
    void continuation(Args&&... args) {
        add_to_target(std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void declare(true_type, Target& target, string name, Args&&... args) {
        continuation(target, name, std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void declare(false_type, Target& target, string, Args&&... args) {
        continuation(target, std::forward<Args>(args)...);
    }
}  // namespace needs_name

// forwarding based on target traits
namespace is_struct {
    // clang-format off
    struct tag {};
    template <class T> struct trait :
        public integral_constant<bool, is_base_of<tag, T>::value> {};
    // clang-format on
}  // namespace is_struct

// forwarding based on target traits
namespace unpack_structs {
    // clang-format off
    struct tag {};
    template <class T> struct trait :
        public integral_constant<bool, is_base_of<tag, T>::value> {};
    // clang-format on

    template <class Target, class... Args>
    void continuation(Target& target, Args&&... args) {
        needs_name::declare(needs_name::trait<Target>(), target, std::forward<Args>(args)...);
    }

    template <class Target, class Variable>
    void declare_helper(true_type, true_type, Target& target, string name, Variable& var) {
        var.struct_interface(target);
    }

    template <class Anything, class... Args>
    void declare_helper(false_type, Anything, Args&&... args) {
        continuation(std::forward<Args>(args)...);
    }

    template <class... Args>
    void declare_helper(true_type, false_type, Args&&... args) {
        continuation(std::forward<Args>(args)...);
    }

    template <class TraitValue, class Target, class Variable, class... Args>
    void declare(TraitValue trait, Target& target, string name, Variable& var, Args&&... args) {
        declare_helper(
            trait, is_struct::trait<Variable>(), target, name, var, std::forward<Args>(args)...);
    }
}  // namespace unpack_structs


template <class Target, class... Args>
void declare(Target& target, string name, Args&&... args) {
    unpack_structs::declare(
        unpack_structs::trait<Target>(), target, name, std::forward<Args>(args)...);
}


// =========================================================================
//   TOOLS
// =========================================================================

struct MyPrinter : needs_name::tag, unpack_structs::tag {
    void add(string name, double value) {
        p->message("Variable %s has value %f", name.c_str(), value);
    }

    void add(string name, int value) {
        p->message("Variable %s has value %d", name.c_str(), value);
    }
};

// clang-format off
// template <> struct needs_name::trait<MyPrinter> : public true_type {};
// clang-format on


// =========================================================================
//   MODEL
// =========================================================================

struct MyStruct : is_struct::tag {
    int j = 3;
    double t = 9.2;

    template <class Target>
    void struct_interface(Target& target) {
        declare(target, "j", j);
        declare(target, "t", t);
    }
};

class MyModel {
    double x = 2.3;
    int i = 5;
    MyStruct s;

  public:
    template <class Target>
    void interface(Target& target) {
        declare(target, "x", x);
        declare(target, "i", i);
        declare(target, "s", s);
    }
};


void compute(int, char**) {
    MyModel m;
    MyPrinter p;
    m.interface(p);
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }