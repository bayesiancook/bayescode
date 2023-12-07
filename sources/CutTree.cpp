#include "Tree.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char* argv[])    {

    Tree tree(argv[1]);
    double cutoff = atof(argv[2]);
    string basename = argv[3];
    tree.SetIndices();
    ofstream os((basename + ".taxgroups").c_str());
    tree.GetRecentGroups(cutoff, os);
}

