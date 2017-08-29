
#include "Random.hpp"
#include "IncompleteGamma.hpp"

int main(int argc, char* argv[])    {

    if (argc == 1)  {
        cerr << "chi2 z df\n";
        cerr << "\tdf = 1 or 2: number of degrees of freedom for simple chi2\n";
        cerr << "\tdf = -1: mixture of 50% point mass 50% chi2 df = 1\n";
        cerr << "\tdf = -2: mixture of 25% point mass, 50% chi2 df=1, 25% chi2 df=2\n";
        exit(1);
    }

    // 1,2 : 1 or 2 degrees of freedom
    // -1 : mix of 50% 0 50% 1 
    // -2: mix of 25% 0 50% 1 25% 2
    double x = atof(argv[1]);
    int df = atoi(argv[2]);

    cout << Chi2Pval(x,df) << '\n';
}

