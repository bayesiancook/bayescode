#include "tree-implem.hpp"
#include <fstream>

int main() {
    std::ifstream file("../../nhx-parser/data/tree1.nhx");
    NHXParser parser{file};
    auto tree = make_from_parser(parser);
    std::cout << "Nb nodes = " << tree->nb_nodes() << std::endl;
}
