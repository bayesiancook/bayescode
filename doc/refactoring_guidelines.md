# New style for suffstats

Compared to legacy bayescode, there is a new way to package suffstats.

Advantages of this new style are:
* **(important)** in debug mode, there is a runtime check that the value is the same before and after gathering;
* homogeneous interface for all suffstats.

Here is the gist of it:
* If a suffstat provides several summary values (e.g. `count` and `beta` for omega path suffstats) they should declare a struct with all summary values.
  * this struct should provide an equality operator (i.e. `operator==`);
  * this struct is meant to be copied; if summary values are expensive to copy then the struct should hold references to the values instead.
* The suffstat object must inherit from `Proxy<T>` where T is the summary value struct (or just the value type if there is a single value).
  * member functions `T _get()` and `void gather()` must be provided;
  * both should have the `final` keyword;
  * the suffstat object itself should be final;
  * `T _get()` should be private.

Here is an example for an omega suffstat wrapping legacy bayescode objects:

```cpp
#include "submodels/suffatt_wrappers.hpp"

struct omega_suffstat_t {
    int count;
    double beta;
    bool operator==(const omega_suffstat_t& other) const {
        return count == other.count && beta == other.beta;
    }
};

class OmegaSSW final : public Proxy<omega_suffstat_t> {  // SSW = suff stat wrapper
    const OmegaCodonSubMatrix& _codon_submatrix;
    const PathSuffStat& _path_suffstat;
    OmegaPathSuffStat _ss;

    omega_suffstat_t _get() final { return {_ss.GetCount(), _ss.GetBeta()}; }

  public:
    OmegaSSW(const OmegaCodonSubMatrix& codon_submatrix, const PathSuffStat& pathsuffstat)
        : _codon_submatrix(codon_submatrix), _path_suffstat(pathsuffstat) {}

    void gather() final {
        _ss.Clear();
        _ss.AddSuffStat(_codon_submatrix, _path_suffstat);
    }
};
```

Functions which need the suffstat value(s) should work with a `Proxy<T>&`. Here is an example use of a suffstat:

```cpp
template <class GlobomModel, class Gen>
static void gibbs_resample( // using suffstat interface reference
    GlobomModel& model, Proxy<omega_suffstat_t>& ss, Gen& gen) {
    double alpha = get<omega, params, shape>(model)();
    double beta = 1. / get<omega, params, struct scale>(model)();
    auto ss_value = ss.get(); // get sufftat value (local copy)
    get<omega, value>(model) = // access ss value fields
        gamma_sr::draw(alpha + ss_value.count, beta + ss_value.beta, gen);
}
```