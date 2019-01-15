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


template <class _Context, class... Args>
void declare(_Context& context, Args&&... args) {
    ignore_name::declare(typename has_tag<typename _Context::context, ignore_name::tag>::type(),
        context, std::forward<Args>(args)...);
}

template <class... Tags, class Target, class Source>
void use_interface(Target& target, Source& source) {
    source.interface(make_decl_info<Tags...>(target));
}


class MyModel {
    int i = 3;
    double j = 2.3;

  public:
    template <class _Context>  // TODO: change name so that Context can be used here
    void interface(_Context context) {
        declare(context, "i", i);
        declare(context, "j", j);
    }
};

struct User {
    template <class T>
    void add(string name, T& i) {
        cout << "Variable " << name << " has value " << i << endl;
    }
};

struct User2 {
    template <class T>
    void add(T& i) {
        cout << "Variable has value " << i << endl;
    }
};

int main() {
    User u;
    User2 u2;
    MyModel m;
    use_interface(u, m);
    use_interface<ignore_name::tag>(u2, m);
}