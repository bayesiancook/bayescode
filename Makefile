# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: _build
	@cd _build ; make --no-print-directory -j8

_build: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

_build_coverage: CMakeLists.txt
	@rm -rf _build_coverage
	@mkdir _build_coverage
	@cd _build_coverage ; cmake -DCOVERAGE_MODE=ON ..

.PHONY: coverage
coverage: _build_coverage
	@make --no-print-directory test

.PHONY: clean
clean:
	@rm -rf _build
	@rm -rf _build_coverage
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
POINTS=5

.PHONY: run-unit-tests
run-unit-tests: all
	@echo "\n\e[35m\e[1m== Tree test ==================================================================\e[0m"
	_build/tree_test
	@echo "\n\n\e[35m\e[1m== All sequential tests =======================================================\e[0m"
	_build/all_tests
	@echo "\n\n\e[35m\e[1m== MPI par test ===============================================================\e[0m"
	mpirun -np 3 _build/mpi_par_test

.PHONY: run-app-tests
run-app-tests: all
	@rm -rf _test
	@mkdir _test
	@echo "\n\e[35m\e[1m== Globom run ===============================================================\e[0m"
	_build/globom -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom restart ===========================================================\e[0m"
	_build/globom _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom read ==============================================================\e[0m"
	_build/readglobom _test/globom_gal4
	@echo "\n\e[35m\e[1m== CodonM2a run =============================================================\e[0m"
	_build/codonm2a -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a restart =========================================================\e[0m"
	_build/codonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a read ============================================================\e[0m"
	_build/readcodonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== MutSel run ===============================================================\e[0m"
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel restart ===========================================================\e[0m"
	_build/aamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read ==============================================================\e[0m"
	_build/readaamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read site-profiles ================================================\e[0m"
	_build/readaamutsel --ss _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== Multigene Single Omega ===================================================\e[0m"
	cd data/small_multigene && mpirun -np 2 ../../_build/multigeneglobom -t tree.nwk -a verysmall.list  -u 1 tmp

.PHONY: test
test:
	@make --no-print-directory run-unit-tests
	@make --no-print-directory run-app-tests

.PHONY: aamutsel
aamutsel: _build
	@cd _build ; make --no-print-directory -j8 aamutsel readaamutsel
	@rm -rf _aamutsel
	@mkdir _aamutsel
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 _aamutsel/gal4
	_build/aamutsel _aamutsel/gal4
	_build/readaamutsel _aamutsel/gal4
	_build/readaamutsel --ss _aamutsel/gal4
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 -p _aamutsel/gal4_poly
	_build/aamutsel _aamutsel/gal4_poly
	_build/readaamutsel _aamutsel/gal4_poly
	_build/readaamutsel --ss _aamutsel/gal4_poly
