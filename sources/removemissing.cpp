
#include <fstream>
#include "CodonSequenceAlignment.hpp"
#include "Random.hpp"

int main(int argc, char *argv[]) {
    string datafile = argv[1];
    SequenceAlignment *data = new FileSequenceAlignment(datafile);
    CodonSequenceAlignment *codondata = new CodonSequenceAlignment(data, Universal);

    ofstream os(argv[2]);
    codondata->ToStreamTriplet(os);
}
