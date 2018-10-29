#include <fstream>
#include <map>
#include "implem.hpp"

int tree_size(const Tree* tree, Tree::NodeIndex from) {
    if (tree->is_leaf(from)) { return 1; }
    int tot = 0;
    for (auto c : tree->children(from)) { tot += tree_size(tree, c); }
    return tot;
}

int main() {
    // parsing file
    std::ifstream file("data/besnard/cyp_coding.Chrysithr_root.nhx");
    NHXParser parser{file};

    // creating topology
    auto tree = make_from_parser(parser);
    std::cout << "Nb nodes = " << tree->nb_nodes() << std::endl;
    std::cout << "root index = " << tree->root() << '\n';

    const Tree* tree2 = tree.get();

    std::cout << "tree size : " << tree_size(tree2, tree2->root()) << '\n';

    for (size_t i = 0; i < tree->nb_nodes(); i++) {
        std::cout << "[" << i << ": " << tree->node_name(i) << "] ";
    }
    std::cout << "\n\n";

    // creating vector of node names
    auto node_names = node_container_from_parser<std::string>(
        parser, [](int i, const AnnotatedTree& t) { return t.tag(i, "name"); });
    for (auto e : node_names) {
        if (e != "") std::cout << e << " ";
    }
    std::cout << "\n\n";

    // slightly shuffled taxa list for example tree
    std::vector<std::string> taxa = {"Ele.bald", "Ele.bal2", "Ele.bal4", "Ele.vivi", "Ele.vivA",
        "Ele.bal3", "Ele.viv2", "Ele.fici", "Ele.grac", "Ele.lim2", "Ele.rost", "Ele.limo",
        "Ele.pal2", "Ele.acut", "Ele.palu", "Ele.gra2", "Ele.lim3", "Ele.geni", "Ele.quan",
        "Abildgaar", "Bulbostyl", "Bolboscho", "Fuir.abn", "Fuir.umb", "Scho.lac", "Scho.val",
        "Scho.muc", "Hellmut1", "Actinosch", "Fimb.lit", "Fimb.dic", "Fimb.fe2", "Fimb.li2",
        "Fimb.di2", "Fimb.fer", "Isolepis", "Hellmut2", "Scirpoid", "Cyp.spha", "Cyp.alt3",
        "Cyp.era6", "Cyp.era1", "Eriophor", "Scirpus", "Schoenox", "Uncin.un", "Uncin.ph",
        "Carex.com", "Carex.hal", "Cyp.fusc", "Cyp.pulc", "Cyp.capi", "Volkiell", "Cyp.ust2",
        "Remirea", "Cyp.iria", "Killinga", "Pycreus", "Cyp.long", "Cyp.rotu", "Cyp.papy",
        "Cyp.ustu", "Blysmus", "Carex.ber", "Carex.pen", "Rhy.alba", "Rhy.grac", "Rhy.albi",
        "Rhy.rubr", "Rhy.glob", "Rhy.glo2", "Carpha", "Schoenus", "Baumea", "Machaeri", "Cladium",
        "Coleochlo", "Microdra", "Chrysithr"};
    auto taxa_index = taxa_index_from_parser(parser, taxa);
    for (auto i : taxa_index) { std::cout << i << " "; }
    std::cout << "\n\n";

    auto taxa_conditions = taxa_container_from_parser<int>(
        parser, taxa, [](int i, const AnnotatedTree& t) { return stoi(t.tag(i, "Condition")); });
    for (size_t i = 0; i < taxa_conditions.vector().size(); i++) {
        std::cout << "[" << taxa.at(i) << ", cond:" << taxa_conditions.vector().at(i)
                  << ", top:" << taxa_conditions.topology_index(i) << "] ";
    }
    std::cout << "\n";
}
