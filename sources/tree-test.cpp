#include <fstream>
#include "tree-implem.hpp"

int main() {
    // parsing file
    std::ifstream file("../data/cyp_coding.nhx");
    NHXParser parser{file};

    // creating topology
    auto tree = make_from_parser(parser);
    std::cout << "Nb nodes = " << tree->nb_nodes() << std::endl;

    // creating vector of node names
    auto node_names = node_container_from_parser<std::string>(
        parser, [](int i, const AnnotatedTree& t) { return t.tag(i, "name"); });
    for (auto e : node_names) {
        if (e != "") std::cout << "Node name: " << e << std::endl;
    }

    // slightly shuffled taxa list for example tree
    std::vector<std::string> taxa = {
        "Ele.bald",  "Ele.bal2",  "Ele.bal4", "Ele.vivi", "Ele.vivA", "Ele.bal3",  "Ele.viv2",
        "Ele.fici",  "Ele.grac",  "Ele.lim2", "Ele.rost", "Ele.limo", "Ele.pal2",  "Ele.acut",
        "Ele.palu",  "Ele.gra2",  "Ele.lim3", "Ele.geni", "Ele.quan", "Abildgaar", "Bulbostyl",
        "Bolboscho", "Fuir.abn",  "Fuir.umb", "Scho.lac", "Scho.val", "Scho.muc",  "Hellmut1",
        "Actinosch", "Fimb.lit",  "Fimb.dic", "Fimb.fe2", "Fimb.li2", "Fimb.di2",  "Fimb.fer",
        "Isolepis",  "Hellmut2",  "Scirpoid", "Cyp.spha", "Cyp.alt3", "Cyp.era6",  "Cyp.era1",
        "Eriophor",  "Scirpus",   "Schoenox", "Uncin.un", "Uncin.ph", "Carex.com", "Carex.hal",
        "Cyp.fusc",  "Cyp.pulc",  "Cyp.capi", "Volkiell", "Cyp.ust2", "Remirea",   "Cyp.iria",
        "Killinga",  "Pycreus",   "Cyp.long", "Cyp.rotu", "Cyp.papy", "Cyp.ustu",  "Blysmus",
        "Carex.ber", "Carex.pen", "Rhy.alba", "Rhy.grac", "Rhy.albi", "Rhy.rubr",  "Rhy.glob",
        "Rhy.glo2",  "Carpha",    "Schoenus", "Baumea",   "Machaeri", "Cladium",   "Coleochlo",
        "Microdra",  "Chrysithr"};
    auto taxa_index = taxa_index_from_parser(parser, taxa);
    for (auto i : taxa_index) {
        std::cout << i << " ";
    }
    std::cout << "\n";
}