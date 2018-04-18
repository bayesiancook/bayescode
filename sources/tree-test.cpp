#include "tree-implem.hpp"

int main() {
    auto tree = make_from_nhx_file("../../nhx-parser/data/tree1.nhx");
    std::cout << "Nb nodes = " << tree->nb_nodes() << std::endl;
}
