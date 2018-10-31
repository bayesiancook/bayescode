# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
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

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================
.PHONY: test
test: cmake
	@cd _build ; make --no-print-directory -j8 tree_test mpi_seq_test mpi_par_test
	@echo "\n\e[35m\e[1m== Tree test ==================================================================\e[0m"
	_build/tree_test
	@echo "\n\n\e[35m\e[1m== MPI seq test ===============================================================\e[0m"
	_build/mpi_seq_test
	@echo "\n\n\e[35m\e[1m== MPI par test ===============================================================\e[0m"
	mpirun -np 3 _build/mpi_par_test

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
