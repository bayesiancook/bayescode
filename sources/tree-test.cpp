#include <fstream>
#include "tree-implem.hpp"

int main() {
    // parsing file
    std::ifstream file("../../nhx-parser/data/tree1.nhx");
    NHXParser parser{file};

    // creating topology
    auto tree = make_from_parser(parser);
    std::cout << "Nb nodes = " << tree->nb_nodes() << std::endl;

    // creating vector of node names
    auto node_names = node_container_from_parser<std::string>(
        parser, [](int i, const AnnotatedTree& t) { return t.tag(i, "name"); });
    for (auto e : node_names) {
        std::cout << "Node name: " << e << std::endl;
    }
}
