# BayesCode

[![Build Status](https://travis-ci.org/bayesiancook/bayescode.svg?branch=dev)](https://travis-ci.org/bayesiancook/bayescode)
[![codecov](https://codecov.io/gh/bayesiancook/bayescode/branch/dev/graph/badge.svg)](https://codecov.io/gh/bayesiancook/bayescode)

If you do not want to compile _BayesCode_, the precompiled binaries for `nodemutsel`, `readnodemutsel`, `mutselomega`
and `readmutselomega` are available [here](https://github.com/ThibaultLatrille/bayescode/releases).
You can then skip the next section.

## How to download and build

To get _BayesCode_ from a machine connected to the internet, type in a terminal:
```bash
git clone https://github.com/ThibaultLatrille/bayescode
```

This should create a folder called `bayescode` (the _BayesCode_ root folder). You must go there before compiling _BayesCode_:
```bash
cd bayescode
```

Requirements to compile _BayesCode_: Clang (or g++) and compatible MPI compiler
```bash
sudo apt install make cmake clang openmpi-bin openmpi-common libopenmpi-dev
```

Then, to build _BayesCode_ simply run:
```bash
make release
```

To check that everything ran well, look into the `bin`  folder to see if executables are present:
```bash
ls bin
```

This should display a list of files. If you see the following files, then _BayesCode_ was correctly built:
* `nodemutsel`
* `readnodemutsel`
* `mutselomega`
* `readmutselomega`

To update to the latest version of _BayesCode_ simply run:
```bash
git pull
make clean
make release
```

## How to format your data

_BayesCode_ requires a tree and an alignment file to run.

Your alignment file must follow the Phylip format.
In addition, the number of bases should be a multiple of 3 (as it will be interpreted as codons).
For example, the following file would be a valid alignment for _BayesCode_ (`nodemutsel` or `mutselomega`):

```phylip
8 6
S0      TCCTGA
S1      AATAGT
S2      GGATTT
S3      AATTCA
S4      CGAAGG
S5      AACGCT
S6      ACGAGT
S7      AATATT
```
A python3 script (no dependencies) to convert Fasta to Phylip is available in the folder `utils`.
```bash
python3 fasta_to_ali.py --input ENSG00000000457_SCYL3_NT.fasta --output ENSG00000000457_SCYL3_NT.phy
```

Your tree file must follow the newick format.
The tree does not need to have branch lenghts.
In addition, the leaves of the tree should have the same names as the sequences in your alignment file.
For example, the following file would be a valid tree file for _BayesCode_ matching the alignment file above :

```newick
((((((((S0,S1),(S2,S3)),(S4,S5),(S6,S7))),(S8,S9),(S10,S11)),(S12,S13),(S14,S15))))
```

The `datà` folder in the _BayesCode_ root folder contains examples of data files usable with _BayesCode_.

## How to run BayesCode

To get the help of the program and the possible options, type:
```bash
bin/nodemutsel --help
```

Basic usage for _BayesCode_ is (from the _BayesCode_ root folder):
```bash
bin/nodemutsel -t <path to tree file> -a <path to alignment file> [options...] <name of run>
```
The name of the run is a string which is used to create files to store the results of the run.
Note that the name of the run must be at the very end of the command.

To restart a chain that was aborted, you just need to give the name of the chain:
```bash
bin/nodemutsel <name of run>
```

## Site-specific adaptive evolution

If you want to use the Mutation-Selection framework for detecting site-specific adaptive evolution, use the program `mutselomega` to run the MCMC and `readmutselomega` to read the trace.

**I. Mutation-selection codon models**

Mutation-selection codon models are obtained by running `mutselomega` for 2000 points of MCMC with the options:
```bash
bin/mutselomega --omegashift 0.0 --ncat 30 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --until 2000 run_mutsel_bglobin
```

The gene-specific mutation matrix (μ) is obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) with the options:
```bash
bin/readmutselomega --every 1 --until 2000 --burnin 1000 --nuc run_mutsel_bglobin
```

The site-specific predicted rate of evolution (ω<sub>0</sub><sup>(i)</sup> : a scalar between 0 and 1 for each site) are obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_mutsel_bglobin.ci0.025.tsv` with the options:
```bash
bin/readmutselomega --every 1 --until 2000 --burnin 1000 --confidence_interval 0.025 --omega_0 run_mutsel_bglobin
```
This will output the mean posterior ω<sub>0</sub><sup>(i)</sup> over the MCMC as well as the 95% (1 - 0.025*2) credible interval.
The value of confidence_interval used as input is considered for each side of the distribution, hence a factor 2 for the credible interval.

The collection of site-specific fitness profiles (F<sup>(i)</sup>: 20x1 vector for each site) are then obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_mutsel_bglobin.siteprofiles` with the options:
```bash
bin/readmutselomega --every 1 --until 2000 --burnin 1000 --ss run_mutsel_bglobin
```

A python script `fitness_to_selcoeff.py` (requirements: numpy and pandas) is available in the folder `utils` to compute selection coefficients (S) from a list of codon transitions.\
The script takes as input a collection of site-specific fitness profiles for a gene generated by _BayesCode_, and a list of positions with the associated codon transitions, it will output the list of selection coefficients for each codon transition.
```bash
python3 utils/fitness_to_selcoeff.py --input_transitions data/bglobin/bglobin.transitions.tsv --input_profiles run_mutsel_bglobin.siteprofiles --output bglobin.transitions.selcoeffs.tsv
```
Help for the script is also available:
```bash
python3 utils/fitness_to_selcoeff.py --help
```

**II. Classical (ω-based) codon models**

Classical (ω-based) codon models are obtained by running `mutselomega` for 2000 points of MCMC with the options:
```bash
bin/mutselomega --omegashift 0.0 --freeomega --omegancat 30 --flatfitness -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre -u 2000 run_classical_bglobin
```

The site-specific predicted rate of evolution (ω<sup>(i)</sup> : a positif scalar) is obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_classical_bglobin.ci0.025.tsv` with the options:
```bash
bin/readmutselomega --every 1 --until 2000 --burnin 1000 --confidence_interval 0.025 run_classical_bglobin
```
This will output the mean posterior ω<sup>(i)</sup> over the MCMC as well as the 95% (1 - 0.025*2) credible interval.
The value of confidence_interval used as input is considered for each side of the distribution, hence a factor 2 for the credible interval.

If you use this program, please cite:\
Nicolas Rodrigue, Thibault Latrille, Nicolas Lartillot,\
A Bayesian Mutation–Selection Framework for Detecting Site-Specific Adaptive Evolution in Protein-Coding Genes,\
_Molecular Biology and Evolution_,
Volume 38, Issue 3, March 2021, Pages 1199–1208,\
https://doi.org/10.1093/molbev/msaa265

## Inferring long-term effective population size

If you want to use the Mutation-Selection framework for inferring long-term changes in effective population size (_N<sub>e</sub>_) and mutation rate (_μ_), use the program `nodemutsel` to run the MCMC and `readnodemutsel` to read the trace.

For example, to run a _BayesCode_ chain called `run_nodemutsel_placentalia` on example placental mammals data (from the `data` folder), until 100 points have been written to disc (from the _BayesCode_ root folder):
```bash
bin/nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx -u 100 run_nodemutsel_placentalia
```

**I. Including life-history traits (optional)**

To include life-history traits (optional), use the option `--traits` with the path to the file containing the traits (in log-space) in tsv format:
```bash
bin/nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx --traitsfile data/placentalia/plac.log.lht -u 2000 run_nodemutsel_placentalia
```

The file containing the traits (in log-space) must have the following (tab-delimited) format:

| TaxonName    | Maturity |  Mass | Longevity |
|--------------|:--------:|------:|----------:|
| Bos          |   6.47   | 13.52 |      3.21 |
| Rousettus    |   5.55   |  4.71 |      2.91 |
| Spermophilus |   5.92   |  5.60 |      2.10 |
| Artibeus     |   NaN    |  3.98 |      2.95 |

The columns are:
 - _TaxonName_: the name of the taxon matching the name in the alignment and the tree.
 - As many columns as traits, without spaces or special characters in the trait. 

The values must follow the conditions:
 - be in log-space, meaning they are log transformed, i.e. y=log(x).
 - can be `NaN` to indicate that the trait is not available for that taxon.

A python script `traits_coevol_to_mutsel.py` (requirements: numpy and pandas) is available in the folder `utils` to convert the traits file from _CoEvol_ format (in natural-space, see [github.com/bayesiancook/coevo](https://github.com/bayesiancook/coevol) for the format definition) to the format used by  _BayesCode_ (in log-space):
```bash
python3 utils/traits_coevol_to_mutsel.py --input data/placentalia/plac.lht --output data/placentalia/plac.log.lht
```

**II. Including fossil calibration (optional)**

To include fossil calibrations (optional), use the option `--fossils` with the path to the file containing the calibrations in tsv format:
```bash
bin/nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx --fossils data/placentalia/plac.calibs.tsv -u 2000 run_nodemutsel_placentalia
```

The file containing the calibrations must have the following (tab-delimited) format:

| NodeName         | Age  | LowerBound | UpperBound |
|------------------|:----:|-----------:|-----------:|
| EquusCeratTapir  | 56.0 |         54 |       58.0 |
| SuBoTrHiMeTuLaVi | 65.0 |          0 |       65.0 |
| MusRattus        | 12.0 |         12 |        inf |
| MegadPteroRouse  | 51.5 |         43 |       60.0 |

The columns are:
- _NodeName_: the name of the internal node matching the name in the input tree. `Root` can be used to designate the root of the tree.
- _Age_: the age of the node, must be between _LowerBound_ and _UpperBound_.
- _LowerBound_: the lower bound (youngest) for the fossil calibration (can be 0).
- _UpperBound_: the upper bound (oldest) for the fossil calibration, can be `inf` to indicate that the calibration is not available for that node.

A python script `calibs_coevol_to_mutsel.py` (requirements: ete3 and pandas) is available in the folder `utils` to convert the traits file from the _CoEvol_ format (using the most recent ancestor between two extent taxa, see [github.com/bayesiancook/coevo](https://github.com/bayesiancook/coevol) for the format definition) to the format used by _BayesCode_ (requires a name for internal nodes):
```bash
```bash
python3 utils/calibs_coevol_to_mutsel.py --input data/placentalia/plac.calibs --tree data/placentalia/plac.calibs --output data/placentalia/plac.calibs.tsv
```
If the tree provided does not include name for internal nodes, the script will give them names and re-write the tree in the same file.

**III. Read annotated trees**

To obtain the annotated newick tree (_N<sub>e</sub>_, _μ_, life-history traits if included) from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
bin/readnodemutsel --burnin 1000 --until 2000 --newick run_nodemutsel_placentalia
```
This command will generate a tree in newick extended format (.nhx) for each trait in the natural-space.

A python script `plot_tree.py` (requirements: matplotlib, ete3, numpy and pandas) is available in the folder `utils` to draw the trees (in .pdf) generated by `readnodemutsel`, using the chain name (here `run_nodemutsel_placentalia`):
```bash
python3 utils/plot_tree.py --input run_nodemutsel_placentalia
```

**IV. Read covariance matrix**

To obtain the covariance, correlation and posterior probabilities matrices (_N<sub>e</sub>_, _μ_, life-history traits) from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
bin/readnodemutsel --burnin 1000 --until 2000 --cov run_nodemutsel_placentalia
```
This command will generate the file `run_nodemutsel_placentalia.cov` containing the covariance matrices.

The entries of the matrices are in the order specified in the header.
The posterior probabilities of a positive correlation (pp) are particularly important.
A pp close to 1 means a strong statistical support for a positive correlation, and a pp close to 0 a supported negative correlation (the posterior probabilities for a negative correlation are given by 1 − pp).

**V. Fitness profiles**

To obtain the fitness profiles in the file `run_nodemutsel_placentalia.profiles` from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
bin/readnodemutsel --burnin 1000 --until 2000 --ss run_nodemutsel_placentalia
```

If you use this program, please cite:\
Thibault Latrille, Vincent Lanore, Nicolas Lartillot,\
Inferring Long-Term Effective Population Size with Mutation–Selection Models,\
_Molecular Biology and Evolution_,
Volume 38, Issue 10, October 2021, Pages 4573–4587,\
https://doi.org/10.1093/molbev/msab160

Moreover, the repository at https://github.com/ThibaultLatrille/MutationSelectionDrift is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript (please use [_BayesCode v1.0_](https://github.com/ThibaultLatrille/bayescode/releases/tag/v1.0)).

### Authors

Nicolas Lartillot (https://github.com/bayesiancook), Thibault Latrille (https://github.com/ThibaultLatrille), Vincent Lanore (https://github.com/vlanore) and Philippe Veber (https://github.com/pveber)
