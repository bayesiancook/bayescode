# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: cmake
	@cd _build ; make --no-print-directory -j8

.PHONY: cmake
cmake: build_dir
	@cd _build ; cmake ..

build_dir: CMakeLists.txt
	@rm -rf _build
	@mkdir _build

.PHONY: clean
clean:
	@rm -rf _build
	@rm -rf _test
	@rm -rf _aamutsel

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================
POINTS=20

.PHONY: run-test
run-test: all
	@rm -rf _test
	@mkdir _test
	@echo "\n\e[35m\e[1m== Tree test ================================================================\e[0m"
	_build/tree_test
	@echo "\n\e[35m\e[1m== MPI seq test =============================================================\e[0m"
	_build/mpi_seq_test
	@echo "\n\e[35m\e[1m== MPI par test =============================================================\e[0m"
	mpirun -np 3 _build/mpi_par_test
	@echo "\n\e[35m\e[1m== Globom run ===============================================================\e[0m"
	_build/globom -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom restart ===========================================================\e[0m"
	_build/globom _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom read ==============================================================\e[0m"
	_build/readglobom -b 0 -e 1 -u ${POINTS} _test/globom_gal4
	@echo "\n\e[35m\e[1m== CodonM2a run =============================================================\e[0m"
	_build/codonm2a -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a restart =========================================================\e[0m"
	_build/codonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a read ============================================================\e[0m"
	_build/readcodonm2a -b 0 -e 1 -u ${POINTS} _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== MutSel run ===============================================================\e[0m"
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel restart ===========================================================\e[0m"
	_build/aamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read ==============================================================\e[0m"
	_build/readaamutsel -b 0 -e 1 -u ${POINTS} _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read site-profiles ================================================\e[0m"
	_build/readaamutsel --ss -b 0 -e 1 -u ${POINTS} _test/aamutsel_gal4

.PHONY: test
test: cmake
	@make run-test

.PHONY: coverage
coverage: build_dir
	@cd _build ; cmake -DCOVERAGE_MODE=ON ..
	@make run-test

.PHONY: aamutsel
aamutsel: cmake
	@cd _build ; make --no-print-directory -j8 aamutsel readaamutsel
	@rm -rf _aamutsel
	@mkdir _aamutsel
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 _aamutsel/gal4
	_build/aamutsel _aamutsel/gal4
	_build/readaamutsel -b 0 -e 1 -u 10 _aamutsel/gal4
	_build/readaamutsel --ss -b 0 -e 1 -u 10 _aamutsel/gal4
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 -p _aamutsel/gal4_poly
	_build/aamutsel _aamutsel/gal4_poly
	_build/readaamutsel -b 0 -e 1 -u 10 _aamutsel/gal4_poly
	_build/readaamutsel --ss -b 0 -e 1 -u 10 _aamutsel/gal4_poly
