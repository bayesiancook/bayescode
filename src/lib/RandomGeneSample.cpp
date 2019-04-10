#include "global/Random.hpp"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    ifstream is(argv[1]);
    int P = atoi(argv[2]);
    ofstream os(argv[3]);
    int N;
    is >> N;
    vector<string> gene(N);
    for (int i = 0; i < N; i++) { is >> gene[i]; }

    int *sample = new int[P];
    Random::DrawFromUrn(sample, P, N);
    os << P << '\n';
    for (int i = 0; i < P; i++) { os << gene[sample[i]] << '\n'; }
}
