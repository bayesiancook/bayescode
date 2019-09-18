#pragma once

#include "structure/type_tag.hpp"
#include "tags/traits.hpp"

template <class MetaData, class... Fields>
struct external_interface<tagged_tuple<MetaData, Fields...>> {
    using TTuple = tagged_tuple<MetaData, Fields...>;

    template <class Info, class Field, class Key>
    static void declare_field(node_tag, Info info, Field& field, Key) {
        DEBUG("Adding node {}.", Key::to_string());
        model_node(info, Key::to_string(), get<value>(field));
    }

    template <class Info, class Field, class Key>
    static void declare_field(model_tag, Info info, Field& field, Key) {
        DEBUG("Adding model {}.", Key::to_string());
        model_node(info, Key::to_string(), field);
    }

    template <class Info, class Field, class Key>
    static void declare_field(unknown_tag, Info info, Field& field, Key) {}

    template <class Info, class... Keys>
    static void declare_interface_helper(Info info, TTuple& target, std::tuple<Keys...>) {
        std::vector<int> ignore = {
            (declare_field(type_tag(get<Keys>(target)), info, get<Keys>(target), Keys{}), 0)...};
    }

    template <class Info>
    static void declare_interface(Info info, TTuple& target) {
        declare_interface_helper(info, target, map_key_list_t<field_map_t<TTuple>>{});
    }
};