# BayesCode

Bayescode is an ongoing project with several ramifications.

### Site-specific adaptive evolution
If you want to use the Mutation-Selection framework for detecting site-specific adaptive evolution, please use the following fork of BayesCode: https://github.com/ThibaultLatrille/bayescode. 


- Nicolas Rodrigue, Thibault Latrille, Nicolas Lartillot, A Bayesian Mutation–Selection Framework for Detecting Site-Specific Adaptive Evolution in Protein-Coding Genes, _Molecular Biology and Evolution_, Volume 38, Issue 3, March 2021, Pages 1199–1208, https://doi.org/10.1093/molbev/msaa265


### Inferring long-term effective population size
If you want to use the Mutation-Selection framework for inferring long-term effective population size, please use the following fork of BayesCode: https://github.com/ThibaultLatrille/bayescode. 

Moreover, the repository at https://github.com/ThibaultLatrille/MutationSelectionDrift is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript and gives the tools to produce your own experiment on your dataset, given you have at least an DNA alignment file and an associated rooted tree topology (branch lengths are not required).

- Thibault Latrille, Vincent Lanore, Nicolas Lartillot, Inferring Long-Term Effective Population Size with Mutation–Selection Models, _Molecular Biology and Evolution_, Volume 38, Issue 10, October 2021, Pages 4573–4587, https://doi.org/10.1093/molbev/msab160


## Get BayesCode up and running on Linux

### Clone the repository

```
$ git clone https://github.com/bayesiancook/bayescode.git
```

### Compile BayesCode
Requirements: mpic++
```
$ cd ./bayescode/sources
$ make
```

### Run BayesCode

```
$ cd ../data
$ mpirun -n 4 multigeneaamutselddp
```

### Authors

Nicolas Lartillot (https://github.com/bayesiancook)
