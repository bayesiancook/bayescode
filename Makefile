# COMPILATION
# Requires: cmake 3.1.0 or better
.PHONY: all
all: cmake
	@cd _build ; make --no-print-directory -j8

.PHONY: cmake
cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

.PHONY: clean
clean:
	@rm -rf _build

.PHONY: test
test: cmake
	@cd _build ; make --no-print-directory -j8 tree_test
	@cd _build ; make --no-print-directory -j8 mpi_seq_test
	@cd _build ; make --no-print-directory -j8 mpi_par_test
	@echo "\n== Tree test ====================================================="
	_build/tree_test
	@echo "\n\n== MPI seq test =================================================="
	_build/mpi_seq_test
	@echo "\n\n== MPI par test =================================================="
	_build/mpi_par_test

# CODE QUALITY
# Requires: clang-format
.PHONY: format
format:
	@clang-format -i `find -name *.*pp`

# Run AaMutSel
.PHONY: aamutsel
aamutsel: all
	@rm -f gal4*.*
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 gal4
	_build/aamutsel gal4
	_build/readaamutsel --om -b 0 -e 1 -u 10 gal4
	_build/readaamutsel --ss -b 0 -e 1 -u 10 gal4
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 -p gal4_poly
	_build/aamutsel gal4_poly
	_build/readaamutsel --om -b 0 -e 1 -u 10 gal4_poly
	_build/readaamutsel --ss -b 0 -e 1 -u 10 gal4_poly
