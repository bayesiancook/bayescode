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

    // slightly shuffled taxa list for example tree
    std::vector<std::string> taxa = {"ENSDNOP00000000726",
                                     "ENSLAFP00000012942",
                                     "ENSMPUP00000016079",
                                     "ENSFCAP00000014291",
                                     "ENSMLUP00000004199",
                                     "ENSSSCP00000026727",
                                     "ENSSSCP00000008253",
                                     "APAAD00000001122_ENSGT00390000001466",
                                     "ENSOARP00000011537",
                                     "ENSP00000344087",
                                     "ENSCJAP00000044683",
                                     "ENSOGAP00000003817",
                                     "ENSCPOP00000005726",
                                     "APAMA00000001586_ENSGT00390000001466",
                                     "TRAMM00000000589_ENSGT00390000001466",
                                     "APAMU00000001207_ENSGT00390000001466",
                                     "ENSOCUP00000007580",
                                     "ENSDNOP00000032016",
                                     "ENSLAFP00000001205",
                                     "ENSFCAP00000021722",
                                     "ENSOCUP00000012223",
                                     "ENSMPUP00000012148",
                                     "ENSMLUP00000008482",
                                     "ENSOARP00000007413",
                                     "ENSSSCP00000019720",
                                     "ENSOGAP00000010184",
                                     "ENSCJAP00000005006",
                                     "ENSP00000348489",
                                     "ENSCPOP00000017070",
                                     "ENSSTOP00000022694",
                                     "APAMA00000001608_ENSGT00390000001466",
                                     "APAAD00000001152_ENSGT00390000001466",
                                     "APAMM00000001564_ENSGT00390000001466",
                                     "APAAD00000001128_ENSGT00390000001466",
                                     "APAMU00000001167_ENSGT00390000001466",
                                     "APAMM00000001562_ENSGT00390000001466",
                                     "ENSDNOP00000023307",
                                     "ENSDNOP00000020805",
                                     "ENSDNOP00000019147",
                                     "ENSDNOP00000029276",
                                     "ENSLAFP00000018072",
                                     "ENSSTOP00000023161",
                                     "APAMM00000001538_ENSGT00390000001466",
                                     "APAAD00000001117_ENSGT00390000001466",
                                     "APAMM00000001535_ENSGT00390000001466",
                                     "APAMA00000001601_ENSGT00390000001466",
                                     "APAMA00000001591_ENSGT00390000001466",
                                     "ENSCPOP00000014781",
                                     "ENSOCUP00000021203",
                                     "ENSP00000396160",
                                     "ENSCJAP00000040067",
                                     "ENSOGAP00000020375",
                                     "ENSMPUP00000002445",
                                     "ENSFCAP00000015322",
                                     "ENSSSCP00000023152",
                                     "ENSOARP00000012163"};
    auto taxa_index = taxa_index_from_parser(parser, taxa);
    for (auto i: taxa_index) {
        std::cout << i << " ";
    }
    std::cout << "\n";
}
