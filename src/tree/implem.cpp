#include "implem.hpp"
#include <map>

std::vector<int> taxa_index_from_parser(TreeParser& parser, const std::vector<std::string>& taxa) {
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    std::map<std::string, int> mapping;  // node name -> node index mapping
    std::set<std::string> taxa_in_tree;  // set of taxa in tree to check consistency
    for (NodeIndex i = 0; i < NodeIndex(tree.nb_nodes()); i++) {
        auto node_name = tree.tag(i, "name");
        if (node_name != "") {
            mapping[node_name] = i;
            taxa_in_tree.insert(node_name);
        }
    }
    if (taxa_in_tree != std::set<std::string>{taxa.begin(), taxa.end()}) {
        std::cerr
            << "Error: taxa list provided to taxa_index_from_parser does not match taxa in tree.\n";
        exit(1);
    }
    std::vector<int> result(taxa.size(), -1);
    for (size_t i = 0; i < taxa.size(); i++) { result.at(i) = mapping.at(taxa.at(i)); }
    return result;
}

std::unique_ptr<const Tree> make_from_parser(TreeParser& parser) {
    return std::unique_ptr<Tree>(new DoubleVectorTree(parser.get_tree()));
}
