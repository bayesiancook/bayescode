# BayesCode

Requirements: Clang (or g++) and compatible MPI compiler
```bash
sudo apt install g++-5 clang-3.6 libopenmpi-dev openmpi-bin
```

## How to download and build

To get BayesCode from a machine connected to the internet, type in a terminal:

```bash
git clone https://github.com/bayesiancook/bayescode.git
```

This should create a folder called `bayescode` (the BayesCode root folder). You must go there before compiling BayesCode:

```bash
cd bayescode
```

Then, to build BayesCode simply run:

```bash
make
```

To check that everything ran well, look into the `_build`  folder to see if executables are present:

```bash
ls _build
```

This should display a list of files. If you see the following files, then BayesCode was correctly built:
* `globom`
* `multigeneglobom`
* `codonm2a`
* `readcodonm2a`
* `aamutsel`
* `readaamutsel`

## How to format your data

BayesCode requires a tree and an alignment file to run.

Your alignment file must follow the Phylip format.
In addition, the number of bases should be a multiple of 3 (as it will be interpreted as codons).
For example, the following file would be a valid alignment for BayesCode (globom, codonm2a or aamutsel):

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

Your tree file must follow the newick format.
The tree does not need to have branch lenghts.
In addition, the leaves of the tree should have the same names as the sequences in your alignment file.
For example, the following file would be a valid tree file for BayesCode matching the alignment file above :

```newick
((((((((S0,S1),(S2,S3)),(S4,S5),(S6,S7))),(S8,S9),(S10,S11)),(S12,S13),(S14,S15))))
```

The `datà` folder in the BayesCode root folder contains examples of data files usable with BayesCode.

## How to run BayesCode

Basic usage for BayesCode is (from the BayesCode root folder):

```
_build/aamutsel -t <path to tree file> -d <path to alignment file> [options...] <name of run>
```
The name of the run is a string which is used to create files to store the results of the run.
Note that the name of the run must be at the very end of the command.

Useful options include:
* `-x <integer e> <integer u>` to tell BayesCode to write to disc every `e` datapoints and to stop after `u` cycles of `e` datapoints.

For example, to run a BayesCode chain called "myrun" on example data (from the `data` folder), until 10k points have been written to disc (from the BayesCode root folder):

```bash
_build/globom -t data/samhd1.tree -d data/samhd1.ali -x 1 10000 myglobomrun
_build/codonm2a -t data/samhd1.tree -d data/samhd1.ali -x 1 10000 mycodonm2arun
_build/aamutsel -t data/samhd1.tree -d data/samhd1.ali -x 1 10000 myaamutselrun
```

To restart a chain that was aborted, you just need to give the name of the chain:
```bash
_build/globom myglobomrun
_build/codonm2a mycodonm2arun
_build/aamutsel myaamutselrun
```

### Authors

Nicolas Lartillot (https://github.com/bayesiancook)
