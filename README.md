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

### Run BayesCode

```
$ cd ../data
$ mpirun -n 4 multigeneaamutselddp
```

### Authors

Nicolas Lartillot (https://github.com/bayesiancook)
