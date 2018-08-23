#include <vector>
#include "nhx-parser.hpp"
#include "tree-interface.hpp"

// a tree with both a vector of parents and a vector of children
class DoubleVectorTree : public TreeTopology {
    using TreeTopology::NodeIndex;
    std::vector<NodeIndex> parent_;
    std::vector<std::set<NodeIndex>> children_;
    NodeIndex root_{0};

  public:
    DoubleVectorTree(const AnnotatedTree& input_tree) {
        root_ = input_tree.root();
        for (std::size_t i = 0; i < input_tree.nb_nodes(); i++) {
            parent_.push_back(input_tree.parent(i));
            children_.emplace_back(input_tree.children(i).begin(), input_tree.children(i).end());
        }
    }

    const std::set<NodeIndex>& children(NodeIndex node) const final { return children_.at(node); }
    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }
    NodeIndex root() const final { return root_; }
    std::size_t nb_nodes() const final { return parent_.size(); }
    bool is_root(NodeIndex i) const final { return i == root_; }
    bool is_leaf(NodeIndex i) const final { return children_.at(i).size() == 0; }
};

std::unique_ptr<const TreeTopology> make_from_parser(TreeParser& parser) {
    return std::unique_ptr<TreeTopology>(new DoubleVectorTree(parser.get_tree()));
}

template <class Element>
std::vector<Element> node_container_from_parser(TreeParser& parser,
                                                Element (*init)(AnnotatedTree::NodeIndex,
                                                                const AnnotatedTree&)) {
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    std::vector<Element> result;
    for (NodeIndex i = 0; i < NodeIndex(tree.nb_nodes()); i++) {
        result.push_back(init(i, tree));
    }
    return result;
}

template <class Element>
std::vector<Element> branch_container_from_parser(TreeParser& parser,
                                                  Element (*init)(AnnotatedTree::NodeIndex,
                                                                  const AnnotatedTree&)) {
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    std::vector<Element> result;
    for (NodeIndex i = 0; i < NodeIndex(tree.nb_nodes()); i++) {
        if (i != tree.root()) {
            result.push_back(init(i, tree));
        } else {
            result.emplace_back();  // default-constructed element for root
        }
    }
    return result;
}

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
    for (size_t i = 0; i < taxa.size(); i++) {
        result.at(i) = mapping.at(taxa.at(i));
    }
    return result;
}

template <class Element>
class TreeElementVector {
    std::vector<Element> v_;
    std::vector<int> index_;

    template <class E>
    friend TreeElementVector<E> taxa_container_from_parser(
        TreeParser& parser, const std::vector<std::string>& taxa,
        E (*init)(AnnotatedTree::NodeIndex, const AnnotatedTree&));

  public:
    using NodeIndex = AnnotatedTree::NodeIndex;

    std::vector<Element>& vector() { return v_; }
    const std::vector<Element>& vector() const { return v_; }
    NodeIndex topology_index(int vector_index) const { return index_.at(vector_index); }
};

template <class Element>
TreeElementVector<Element> taxa_container_from_parser(TreeParser& parser,
                                                      const std::vector<std::string>& taxa,
                                                      Element (*init)(AnnotatedTree::NodeIndex,
                                                                      const AnnotatedTree&)) {
    TreeElementVector<Element> result;
    result.index_ = taxa_index_from_parser(parser, taxa);
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    for (size_t i = 0; i < taxa.size(); i++) {
        result.v_.push_back(init(NodeIndex(result.index_.at(i)), tree));
    }
    return result;
}
