# BayesCode

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
### Compile BayesCode with meson
You can install meson with `brew install meson` (Mac) or `apt install meson` (Linux).

```
$ cd ./bayescode sources
$ meson build -Dprefix=$HOME/Applications/bayescode
$ ninja -C build install

The binaries will end up as e.g. <prefix>/bin/multigeneaamutselddp.

### Run BayesCode

```
$ cd ../data
$ mpirun -n 4 multigeneaamutselddp
```

### Authors

Nicolas Lartillot (https://github.com/bayesiancook)
