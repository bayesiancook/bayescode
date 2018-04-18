#include <fstream>
#include <vector>
#include "tree-implem.hpp"

// a tree with both a vector of parents and a vector of children
class DoubleVectorTree : public TreeTopology {
    using TreeTopology::NodeIndex;
    std::vector<NodeIndex> parent_;
    std::vector<std::set<NodeIndex>> children_;
    NodeIndex root_{0};

  public:
    DoubleVectorTree(const AnnotatedTree &input_tree) {
        root_ = input_tree.root();
        for (std::size_t i = 0; i < input_tree.nb_nodes(); i++) {
            parent_.push_back(input_tree.parent(i));
            children_.emplace_back(input_tree.children(i).begin(), input_tree.children(i).end());
        }
    }

    const std::set<NodeIndex> &children(NodeIndex node) const final { return children_.at(node); }

    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }

    NodeIndex root() const final { return root_; }

    std::size_t nb_nodes() const final { return parent_.size(); }

    bool is_root(NodeIndex i) const final { return i == root_; }

    bool is_leaf(NodeIndex i) const final { return children_.at(i).size() == 0; }
};

std::unique_ptr<TreeTopology> make_from_parser(TreeParser& parser) {
    return std::unique_ptr<TreeTopology>(new DoubleVectorTree(parser.get_tree()));
}
