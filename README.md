# BayesCode
[![Version](https://anaconda.org/bioconda/bayescode/badges/version.svg)](https://anaconda.org/bioconda/bayescode)
[![Release](https://anaconda.org/bioconda/bayescode/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bayescode)
[![Build Status](https://anaconda.org/bioconda/bayescode/badges/platforms.svg)](https://anaconda.org/bioconda/bayescode)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ThibaultLatrille/bayescode/blob/chronogram/License.MD)
[![codecov](https://codecov.io/gh/bayesiancook/bayescode/branch/dev/graph/badge.svg)](https://codecov.io/gh/bayesiancook/bayescode)


The precompiled binaries are available as the `bayescode` conda package from the channel `bioconda`:

```bash
conda install -c conda-forge -c bioconda bayescode
```

**Wiki and documentation at [github.com/ThibaultLatrille/bayescode/wiki](https://github.com/ThibaultLatrille/bayescode/wiki)**

## Wiki table of contents
- [Installation](https://github.com/ThibaultLatrille/bayescode/wiki/1.-installation)
- [Format your data](https://github.com/ThibaultLatrille/bayescode/wiki/2.-format-your-data)
- [Run BayesCode](https://github.com/ThibaultLatrille/bayescode/wiki/3.-run-bayescode)
- [`mutselomega` - Gene and site-specific rates of evolution (ω, ω<sub>0</sub>, ω<sub>A</sub><sup>phy</sup>, ω<sub>∗</sub>) and selection coefficients (S<sub>0</sub>)](https://github.com/ThibaultLatrille/bayescode/wiki/4.-mutselomega)
- [`nodemutsel` - Changes of effective population size (_N<sub>e</sub>_) and mutation rate (_μ_) along the phylogeny](https://github.com/ThibaultLatrille/bayescode/wiki/5.-nodemutsel)
- [`nodeomega` - Changes of ω along the phylogeny](https://github.com/ThibaultLatrille/bayescode/wiki/6.-nodeomega)
- [`nodetraits` - Test of diversifying selection for a quantitative trait](https://github.com/ThibaultLatrille/bayescode/wiki/7.-nodetraits)
- [Citations](https://github.com/ThibaultLatrille/bayescode/wiki/citations)

## References

- **The preprint for use of `mutselomega` to compute S<sub>0</sub> is**:

T. Latrille, J. Joseph, D. A. Hartasánchez, N. Salamin,\
 Mammalian protein-coding genes exhibit widespread beneficial mutations that are not adaptive, \
_bioRxiv_,\
[doi.org/10.1101/2023.05.03.538864](https://doi.org/10.1101/2023.05.03.538864)

_Scripts and data necessary to reproduce figures at [github.com/ThibaultLatrille/SelCoeff](https://github.com/ThibaultLatrille/SelCoeff)._

- **If you use `mutselomega` to compute ω<sub>0</sub> or ω<sub>A</sub><sup>phy</sup> = ω - ω<sub>0</sub>, please cite**:

Thibault Latrille, Nicolas Rodrigue, Nicolas Lartillot,\
Genes and sites under adaptation at the phylogenetic scale also exhibit adaptation at the population-genetic scale,\
_Proceedings of the National Academy of Sciences_,
Volume 120, Issue 11, March 2023, Pages e2214977120,\
[doi.org/10.1073/pnas.2214977120](https://doi.org/10.1073/pnas.2214977120)

_Scripts and data necessary to reproduce figures at [github.com/ThibaultLatrille/AdaptaPop](https://github.com/ThibaultLatrille/AdaptaPop)._

- **If you use `mutselomega` to compute ω<sub>∗</sub>, please cite**:

Nicolas Rodrigue, Thibault Latrille, Nicolas Lartillot,\
A Bayesian Mutation–Selection Framework for Detecting Site-Specific Adaptive Evolution in Protein-Coding Genes,\
_Molecular Biology and Evolution_,
Volume 38, Issue 3, March 2021, Pages 1199–1208,\
[doi.org/10.1093/molbev/msaa265](https://doi.org/10.1093/molbev/msaa265)

- **If you use `nodemutsel` to compute branch/node specific _N<sub>e</sub>_ and _μ_, or if you use `nodeomega` to compute branch/node specific ω and _μ_ please cite**:

Thibault Latrille, Vincent Lanore, Nicolas Lartillot,\
Inferring Long-Term Effective Population Size with Mutation–Selection Models,\
_Molecular Biology and Evolution_,
Volume 38, Issue 10, October 2021, Pages 4573–4587,\
[doi.org/10.1093/molbev/msab160](https://doi.org/10.1093/molbev/msab160)

_Scripts and data necessary to reproduce figures at [github.com/ThibaultLatrille/MutationSelectionDrift](https://github.com/ThibaultLatrille/MutationSelectionDrift) (using [BayesCode v1.0](https://github.com/ThibaultLatrille/bayescode/releases/tag/v1.0))_.

## Authors
Nicolas Lartillot ([github.com/bayesiancook](https://github.com/bayesiancook)), Thibault Latrille ([github.com/ThibaultLatrille](https://github.com/ThibaultLatrille)), Vincent Lanore ([github.com/vlanore](https://github.com/vlanore)), Philippe Veber ([github.com/pveber](https://github.com/pveber)) and Nicolas Rodrigue.
