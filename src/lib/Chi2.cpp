
#include "IncompleteGamma.hpp"
#include "global/Random.hpp"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        cerr << "chi2 infile df outfile\n";
        cerr << "\tdf = 1 or 2: number of degrees of freedom for simple chi2\n";
        cerr << "\tdf = -1: mixture of 50% point mass 50% chi2 df = 1\n";
        cerr << "\tdf = -2: mixture of 25% point mass, 50% chi2 df=1, 25% chi2 "
                "df=2\n";
        exit(1);
    }

    // 1,2 : 1 or 2 degrees of freedom
    // -1 : mix of 50% 0 50% 1
    // -2: mix of 25% 0 50% 1 25% 2
    ifstream is(argv[1]);
    int df = atoi(argv[2]);
    ofstream os(argv[3]);

    int Ngene;
    is >> Ngene;
    for (int i = 0; i < Ngene; i++) {
        string name;
        double x;
        is >> name >> x;
        x *= 2;
        if (x <= 0) {
            os << name << '\t' << x << '\t' << 1 << '\t' << "100" << '\n';
        } else {
            double p = 0.25 * Chi2Pval(x, df);
            double q = (int)(100 * Ngene * p / (i + 1));
            if (q > 100) { q = 100; }
            os << name << '\t' << x << '\t' << p << '\t' << q << '\n';
        }
    }
}
