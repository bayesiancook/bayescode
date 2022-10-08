# BayesCode
[![Build Status](https://travis-ci.org/bayesiancook/bayescode.svg?branch=dev)](https://travis-ci.org/bayesiancook/bayescode)
[![codecov](https://codecov.io/gh/bayesiancook/bayescode/branch/dev/graph/badge.svg)](https://codecov.io/gh/bayesiancook/bayescode)

## Contents
- [Installation](#installation)
  1. Blue pill with conda
  2. Red pill and compile code
- [Format your data](#format-your-data)
  1. Alignment file
  2. Tree file
  3. Example files
- [Run BayesCode](#run-bayescode)
- [Site-specific adaptive evolution](#site-specific-adaptive-evolution)
  1. Mutation-selection codon models (ω<sub>0</sub>)
  2. Classical codon models (ω)
  3. Mutation-selection codon models with a multiplicative factor (ω<sub>∗</sub>)
  4. Options for `mutselomega` and `readmutselomega`
- [Inferring long-term effective population size](#inferring-long-term-effective-population-size)
  1. Including life-history traits
  2. Including fossil calibration
  3. Read annotated trees
  4. Read covariance matrix
  5. Fitness profiles
  6. Options for `nodemutsel` and `readnodemutsel`
- [Citations](#citations)
- [Authors](#authors) 

## Installation
At this step, you may have to chose between either a blue or a red pill :pill:.

With the blue pill, simply install the `bayescode` conda package from the channel `bioconda`.
Or chose the red pill and install the requirements and compile the C++ code for *BayesCode* in your Debian based OS.

Blue pill is preferred is you want minimal conflict with your local system.
Red pill method is preferred if you plan to tinker with the code.
The two pills are mutually not exclusive, no overdose had ever been observed (though no statistical study had been performed).

### 1. Blue pill with conda
The precompiled binaries for `nodemutsel`, `readnodemutsel`, `mutselomega`, `readmutselomega` and python scripts are available as the `bayescode` conda package from the channel `bioconda` ([https://anaconda.org/bioconda/bayescode](https://anaconda.org/bioconda/bayescode)).

```bash
conda install -c conda-forge -c bioconda bayescode
```
You can then go to the section [Format your data](#format-your-data).

### 2. Red pill and compile code
Requirements to compile _BayesCode_ are a C++ compiler (`clang` or `g++`), `make` and `cmake`:
```bash
sudo apt install make cmake clang
```
To get and build _BayesCode_:
```bash
git clone https://github.com/ThibaultLatrille/bayescode
cd bayescode
make tiny
```

To check that everything ran well, look into the `bin`  folder to see if executables are present:
* `nodemutsel`
* `readnodemutsel`
* `mutselomega`
* `readmutselomega`

To update to the latest version of _BayesCode_ simply run:
```bash
git pull
make clean
make tiny
```

**Optional**: to compile the whole _BayesCode_ suite, a compatible MPI compiler is required
```bash
sudo apt install openmpi-bin openmpi-common libopenmpi-dev
make release
```
## Format your data
_BayesCode_ requires a tree and an alignment file to run.

### 1. Alignment file
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
A python3 script to convert Fasta to Phylip is available:
```bash
python3 fasta_to_ali.py --input ENSG00000000457_SCYL3_NT.fasta --output ENSG00000000457_SCYL3_NT.phy
```

### 2. Tree file
Your tree file must follow the newick format.
The tree does not need to have branch lengths.
In addition, the leaves of the tree should have the same names as the sequences in your alignment file.
For example, the following file would be a valid tree file for _BayesCode_ matching the alignment file above:

```newick
((((((((S0,S1),(S2,S3)),(S4,S5),(S6,S7))),(S8,S9),(S10,S11)),(S12,S13),(S14,S15))))
```

### 3. Example files
The `datà` folder in the _BayesCode_ root folder contains examples of data files usable with _BayesCode_.
The whole folder can be downloaded here: [https://github.com/ThibaultLatrille/bayescode/releases/download/v1.1.6/data.zip](https://github.com/ThibaultLatrille/bayescode/releases/download/v1.1.6/data.zip).

## Run BayesCode
To get the help of the program and the possible options, type:
```bash
nodemutsel --help
```

Basic usage for _BayesCode_ is (from the _BayesCode_ root folder):
```bash
nodemutsel -t <path to tree file> -a <path to alignment file> [options...] <name of run>
```
The name of the run is a string which is used to create files to store the results of the run.
Note that the name of the run must be at the very end of the command.

To restart a chain that was aborted, you just need to give the name of the chain:
```bash
nodemutsel <name of run>
```

The global options are:
```
-a <string>,  --alignment <string>
 (required) File path to alignment (PHYLIP format).

-t <string>,  --tree <string>
 (required) File path to the tree (NHX format).

--force
 Overwrite existing output files.

--until <int> (default: -1)
 Maximum number of (saved) iterations (-1 means unlimited).

--every <int> (default: 1)
 Number of MCMC iterations between two saved point in the trace.
```

To read a chain, you need to give the name of the chain to the sister program `read`:
```bash
readnodemutsel --burnin 1000 --until 2000 --every 1 --trace <name of run>
```
The global options for `read` are:
```
--burnin <int>
 Number of MCMC iterations for the burn-in.

--until <int>
 Maximum number of (saved) iterations (-1 means unlimited).

--every <int>
 Number of MCMC iterations between two saved point in the trace.
 
--trace (default: false)
 Recompute the trace.
```
## Site-specific adaptive evolution
If you want to use the Mutation-Selection framework for detecting site-specific adaptive evolution, use the program `mutselomega` to run the MCMC and `readmutselomega` to read the trace.

### I. Mutation-selection codon models (ω<sub>0</sub>)
Mutation-selection codon models are obtained by running `mutselomega` for 2000 points of MCMC with the options:
```bash
mutselomega --ncat 30 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --until 2000 run_mutsel_bglobin
```

The gene-specific mutation matrix (μ) is obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) with the options:
```bash
readmutselomega --every 1 --until 2000 --burnin 1000 --nuc run_mutsel_bglobin
```

The site-specific predicted rate of evolution (ω<sub>0</sub><sup>(i)</sup>: a scalar between 0 and 1 for each site) are obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_mutsel_bglobin.ci0.025.tsv` with the options:
```bash
readmutselomega --every 1 --until 2000 --burnin 1000 --confidence_interval 0.025 --omega_0 run_mutsel_bglobin
```
This will output the mean posterior ω<sub>0</sub><sup>(i)</sup> over the MCMC as well as the 95% (1 - 0.025*2) credible interval.
The value of confidence_interval used as input is considered for each side of the distribution, hence a factor 2 for the credible interval.

The collection of site-specific fitness profiles (F<sup>(i)</sup>: 20x1 vector for each site) are then obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_mutsel_bglobin.siteprofiles` with the options:
```bash
readmutselomega --every 1 --until 2000 --burnin 1000 --ss run_mutsel_bglobin
```

A python script `fitness_to_selcoeff.py` is available to compute selection coefficients (S) from a list of codon transitions.\
The script takes as input a collection of site-specific fitness profiles for a gene generated by _BayesCode_, and a list of positions with the associated codon transitions, it will output the list of selection coefficients for each codon transition.
```bash
fitness_to_selcoeff.py --input_transitions data/bglobin/bglobin.transitions.tsv --input_profiles run_mutsel_bglobin.siteprofiles --output bglobin.transitions.selcoeffs.tsv
```
Help for the script is also available:
```bash
fitness_to_selcoeff.py --help
```

### II. Classical codon models (ω)
Classical codon models (multiplicative factor ω for non-synonymous substitutions) are obtained by running `mutselomega` for 2000 points of MCMC with the options:
```bash
mutselomega --freeomega --omegancat 30 --flatfitness -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre -u 2000 run_classical_bglobin
```

The site-specific predicted rate of evolution (ω<sup>(i)</sup>: a positif scalar) is obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_classical_bglobin.ci0.025.tsv` with the options:
```bash
readmutselomega --every 1 --until 2000 --burnin 1000 --confidence_interval 0.025 run_classical_bglobin
```
This will output the mean posterior ω<sup>(i)</sup> over the MCMC as well as the 95% (1 - 0.025*2) credible interval.
The value of `confidence_interval` used as input is considered for each side of the distribution, hence a factor 2 for the credible interval.

### III. Mutation-selection codon models with a multiplicative factor (ω<sub>∗</sub>)
Mutation-selection codon models with a ω<sub>∗</sub> multiplicative factor (MutSel-M3, see [https://doi.org/10.1093/molbev/msaa265](https://doi.org/10.1093/molbev/msaa265)) are obtained by running `mutselomega` for 2000 points of MCMC with the options:
```bash
mutselomega --freeomega --omegancat 3 --ncat 30 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre -u 2000 run_mutselM3_bglobin
```

The site-specific posterior probabilities in favor of a value greater than 1 (p(ω<sub>∗</sub> > 1)) are obtained by running `readmutselomega`, reading 1000 points of MCMC (first 1000 are considered as burn-in) and written in the file `run_classical_bglobin.omegappgt1.0` as a .tsv file with the options:
```bash
readmutselomega --every 1 --until 2000 --burnin 1000 --omega_threshold 1.0 run_mutselM3_bglobin
```

### IV. Options for `mutselomega` and `readmutselomega`
The options for `mutselomega` are:
```
--ncat <int> (default: 30)
 Number of components for the amino-acid fitness profiles (truncation
 of the stick-breaking process).

--profiles <string> (default: None)
 File path the fitness profiles (tsv or csv), thus considered fixed.
 Each line must contains the fitness of each of the 20 amino-acid, thus
 summing to one. If same number of profiles as the codon alignment,
 site allocations are considered fixed. If smaller than the alignment
 size, site allocations are computed and `ncat` is given by the number
 of profiles in the file.
   
--flatfitness (default: false)
 Fitness profiles are flattened (and `ncat` equals to 1). This option
 is not compatible with the option `profiles`.
    
--freeomega (default: false)
 ω is allowed to vary (default ω is 1.0). Combined with the
 option `flatfitness`, we obtain the classical, ω-based codon model
 (Muse & Gaut). Without the option `flatfitness`, we obtain the
 mutation-selection codon model with a multiplicative factor (ω⁎).

--omegancat <int> (default: 1)
 Number of components for ω (finite mixture).
      
--omegaarray <string> (default: None)
 File path to ω values (one ω per line), thus considered
 fixed. `freeomega` is overridden to false and `omegancat` equals to the
 number of ω in the file.

--omegashift <double> (default: 0.0)
 Additive shift applied to all ω (0.0 for the general case).
```

You can use one of these options with `readmutselomega` to compute a specific posterior statistic:
```
--ss (default: false)
 Computes the mean posterior site-specific amino-acid equilibrium
 frequencies (amino-acid fitness profiles).

--omega_threshold <double> (default: 1.0)
 Threshold to compute the mean posterior probability that ω⁎ (or ω
 if option `flatfitness` is used in `mutselomega`) is greater than a
 given value.

--confidence_interval <string> (default: None)
 Posterior credible interval for ω (per site and at the gene level).

--omega_0 (default: false)
 Posterior credible interval for ω0 predicted at the
 mutation-selection equilibrium from the fitness profiles (instead of
 ω). To use combined with the option `confidence_interval`.

--nuc (default: false)
 Mean posterior nucleotide matrix.
```

## Inferring long-term effective population size
If you want to use the Mutation-Selection framework for inferring long-term changes in effective population size (_N<sub>e</sub>_) and mutation rate (_μ_), use the program `nodemutsel` to run the MCMC and `readnodemutsel` to read the trace.

For example, to run a _BayesCode_ chain called `run_nodemutsel_placentalia` on example placental mammals data (from the `data` folder), until 100 points have been written to disc (from the _BayesCode_ root folder):
```bash
nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx -u 100 run_nodemutsel_placentalia
```

### I. Including life-history traits (optional)
To include life-history traits (optional), use the option `--traits` with the path to the file containing the traits (in log-space) in tsv format:
```bash
nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx --traitsfile data/placentalia/plac.log.lht -u 2000 run_nodemutsel_placentalia
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

A python script `traits_coevol_to_mutsel.py` is available to convert the traits file from _CoEvol_ format (in natural-space, see [github.com/bayesiancook/coevo](https://github.com/bayesiancook/coevol) for the format definition) to the format used by  _BayesCode_ (in log-space):
```bash
traits_coevol_to_mutsel.py --input data/placentalia/plac.lht --output data/placentalia/plac.log.lht
```

### II. Including fossil calibration (optional)
To include fossil calibrations (optional), use the option `--fossils` with the path to the file containing the calibrations in tsv format:
```bash
nodemutsel --ncat 30 -a data/placentalia/plac.ali -t data/placentalia/plac.nhx --fossils data/placentalia/plac.calibs.tsv -u 2000 run_nodemutsel_placentalia
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

A python script `calibs_coevol_to_mutsel.py` is available to convert the traits file from the _CoEvol_ format (using the most recent ancestor between two extent taxa, see [github.com/bayesiancook/coevo](https://github.com/bayesiancook/coevol) for the format definition) to the format used by _BayesCode_ (requires a name for internal nodes):

```bash
calibs_coevol_to_mutsel.py --input data/placentalia/plac.calibs --tree data/placentalia/plac.calibs --output data/placentalia/plac.calibs.tsv
```
If the tree provided does not include name for internal nodes, the script will give them names and re-write the tree in the same file.

### III. Read annotated trees
To obtain the annotated newick tree (_N<sub>e</sub>_, _μ_, life-history traits if included) from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
readnodemutsel --burnin 1000 --until 2000 --newick run_nodemutsel_placentalia
```
This command will generate a tree in newick extended format (.nhx) for each trait in the natural-space.

A python script `plot_tree.py` is available to draw the trees (in .pdf) generated by `readnodemutsel`, using the chain name (here `run_nodemutsel_placentalia`):
```bash
plot_tree.py --input run_nodemutsel_placentalia
```

### IV. Read covariance matrix
To obtain the covariance, correlation and posterior probabilities matrices (_N<sub>e</sub>_, _μ_, life-history traits) from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
readnodemutsel --burnin 1000 --until 2000 --cov run_nodemutsel_placentalia
```
This command will generate the file `run_nodemutsel_placentalia.cov` containing the covariance matrices.

The entries of the matrices are in the order specified in the header.
The posterior probabilities of a positive correlation (pp) are particularly important.
A pp close to 1 means a strong statistical support for a positive correlation, and a pp close to 0 a supported negative correlation (the posterior probabilities for a negative correlation are given by 1 − pp).

### V. Fitness profiles
To obtain the fitness profiles in the file `run_nodemutsel_placentalia.profiles` from the chain `run_nodemutsel_placentalia`, discarding the first 1000 points:
```bash
readnodemutsel --burnin 1000 --until 2000 --ss run_nodemutsel_placentalia
```

### VI. Options for `nodemutsel` and `readnodemutsel`
The options for `nodemutsel` are:
```
--ncat <int> (default: 30)
 Number of components for the amino-acid fitness profiles (truncation
 of the stick-breaking process).
 
--fossils <string> (default: None)
 File path to the fossils calibration in tsv format with columns
 `NodeName`, `Age, `LowerBound` and `UpperBound`.

--profiles <string> (default: None)
 File path the fitness profiles (tsv or csv), thus considered fixed.
 Each line must contains the fitness of each of the 20 amino-acid, thus
 summing to one. If same number of profiles as the codon alignment,
 site allocations are considered fixed. If smaller than the alignment
 size, site allocations are computed and `ncat` is given by the number
 of profiles in the file.

--traitsfile <string> (default: None)
 File path to the life-history trait (in log-space) in tsv format. The
 First column is `TaxonName` (taxon matching the name in the alignment)
 and the next columns are traits.
```

You can use one of these options with `readnodemutsel` to compute a specific posterior statistic:
```
--newick (default: false)
 Computes the mean posterior node-specific entries of the multivariate
 Brownian process. Each entry of the multivariate Brownian process is
 written in a newick extended (.nhx) format file.

--cov (default: false)
 Computes the mean posterior covariance matrix.

--ss (default: false)
 Computes the mean posterior site-specific amino-acid equilibrium
 frequencies (amino-acid fitness profiles).
 
--profiles <string> (default: None)
 Change the profiles filename if desired, otherwise given by
 {chain_name}.siteprofiles as default.
```



## Citations
- If you use `mutselomega` (ω<sub>∗</sub>, ω<sub>0</sub>, ω), please cite:

Nicolas Rodrigue, Thibault Latrille, Nicolas Lartillot,\
A Bayesian Mutation–Selection Framework for Detecting Site-Specific Adaptive Evolution in Protein-Coding Genes,\
_Molecular Biology and Evolution_,
Volume 38, Issue 3, March 2021, Pages 1199–1208,\
https://doi.org/10.1093/molbev/msaa265

- If you use `nodemutsel` (_N<sub>e</sub>_, _μ_), please cite:

Thibault Latrille, Vincent Lanore, Nicolas Lartillot,\
Inferring Long-Term Effective Population Size with Mutation–Selection Models,\
_Molecular Biology and Evolution_,
Volume 38, Issue 10, October 2021, Pages 4573–4587,\
https://doi.org/10.1093/molbev/msab160

The repository at https://github.com/ThibaultLatrille/MutationSelectionDrift is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript (please use [_BayesCode v1.0_](https://github.com/ThibaultLatrille/bayescode/releases/tag/v1.0)).

## Authors
Nicolas Lartillot (https://github.com/bayesiancook), Thibault Latrille (https://github.com/ThibaultLatrille), Vincent Lanore (https://github.com/vlanore), Philippe Veber (https://github.com/pveber) and Nicolas Rodrigue.
