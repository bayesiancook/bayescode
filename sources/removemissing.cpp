
#include "Random.hpp"
#include "CodonSequenceAlignment.hpp"
#include <fstream>

int main(int argc, char* argv[])    {

    string datafile = argv[1];
    SequenceAlignment* data = new FileSequenceAlignment(datafile);
    CodonSequenceAlignment* codondata = new CodonSequenceAlignment(data,Universal);

    ofstream os(argv[2]);
    codondata->ToStreamTriplet(os);
}
