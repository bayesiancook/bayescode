#include <iostream>
#include "decl_info.hpp"

using namespace std;

namespace ignore_name {
    struct tag {};
    template <class T>
    struct trait : public std::integral_constant<bool, std::is_base_of<tag, T>::value> {};

    template <class _DeclInfo, class Ref>
    void declare(true_type, _DeclInfo& info, string, Ref& ref) {
        info.target.add(ref);
    }

    template <class _DeclInfo, class Ref>
    void declare(false_type, _DeclInfo& info, string name, Ref& ref) {
        info.target.add(name, ref);
    }
}  // namespace ignore_name


template <class _DeclInfo, class... Args>
void declare(_DeclInfo& context, Args&&... args) {
    ignore_name::declare(integral_constant < bool,
        _DeclInfo::context::template has_tag<ignore_name::tag>::type::value or
            ignore_name::trait<typename _DeclInfo::target_type>::value > (),
        context, std::forward<Args>(args)...);
}

// template <class _DeclInfo, class... Args>
// void random_var(_DeclInfo& context, Args&&... args) {

// }

template <class... Tags, class Target, class Source>
void use_interface(Target& target, Source& source) {
    source.interface(make_decl_info<Tags...>(target));
}


class MyModel {
    int i = 3, k = 5;
    double j = 2.3;

  public:
    template <class _DeclInfo>
    void interface(_DeclInfo info) {
        declare(info, "i", i);
        declare(info, "j", j);
        declare(info, "k", k);
    }
};

struct Tracer {
    template <class T>
    void add(string name, T& i) {
        cout << "Variable " << name << " has value " << i << endl;
    }
};

struct Tracer2 : ignore_name::tag {
    template <class T>
    void add(T& i) {
        i++;
        cout << "Variable has value+1 " << i << endl;
    }
};

int main() {
    Tracer u;
    Tracer2 u2;
    MyModel m;
    use_interface(u, m);
    use_interface(u2, m);
}