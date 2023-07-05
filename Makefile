# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: bin
	@cd bin ; make --no-print-directory -j8

bin: CMakeLists.txt # default mode is release
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake ..

.PHONY: coverage
coverage:
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake -DCOVERAGE_MODE=ON ..
	@make --no-print-directory test

.PHONY: debug
debug:
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake -DDEBUG_MODE=ON ..
	@make --no-print-directory

.PHONY: release
release:
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake ..
	@make --no-print-directory
	@rm -rf bin/CMakeCache.txt
	@rm -rf bin/cmake_install.cmake
	@rm -rf bin/CMakeFiles
	@rm -rf bin/Makefile
	@rm -rf bin/*.a

.PHONY: tiny
tiny:
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake -DTINY=ON ..
	@make --no-print-directory

.PHONY: clean
clean:
	@rm -rf bin
	@rm -rf bin_coverage
	@rm -rf _test
	@rm -rf _build
	@rm -rf _aamutsel
	@rm -rf _mutseldm5
	@rm -rf _mutselomega
	@rm -rf _dated

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================
POINTS=2

.PHONY: run-unit-tests
run-unit-tests: all
	@echo "\n\e[35m\e[1m== Tree test ==================================================================\e[0m"
	bin/tree_test


.PHONY: run-unit-tests-mpi
run-unit-tests-mpi: all
	@echo "\n\n\e[35m\e[1m== MPI par test ===============================================================\e[0m"
	mpirun -np 3 bin/mpi_par_test

.PHONY: run-app-tests
run-app-tests: all
	@rm -rf _test
	@mkdir _test
	@echo "\n\e[35m\e[1m== Globom run ===============================================================\e[0m"
	@bin/globom -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom restart ===========================================================\e[0m"
	bin/globom _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom read ==============================================================\e[0m"
	bin/readglobom _test/globom_gal4
	@echo "\n\e[35m\e[1m== CodonM2a run =============================================================\e[0m"
	bin/codonm2a -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a restart =========================================================\e[0m"
	bin/codonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a read ============================================================\e[0m"
	bin/readcodonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== MutSel run ===============================================================\e[0m"
	bin/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel restart ===========================================================\e[0m"
	bin/aamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read ==============================================================\e[0m"
	bin/readaamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read site-profiles ================================================\e[0m"
	bin/readaamutsel --ss _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel with polymorphism run =============================================\e[0m"
	bin/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -p -u ${POINTS} _test/aamutsel_gal4_poly
	@echo "\n\e[35m\e[1m== MutSel Multiple omega run ================================================\e[0m"
	bin/mutselomega -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --freeomega --omegancat 3 -u ${POINTS} _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== MutSel Multiple omega restart ============================================\e[0m"
	bin/mutselomega _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== MutSel Multiple omega read ===============================================\e[0m"
	bin/readmutselomega _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== Diffsel double sparse ====================================================\e[0m"
	@make --no-print-directory diffseldsparse
	@echo "\n\e[35m\e[1m== Node Omega run ===========================================================\e[0m"
	bin/nodeomega -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/nodeomega_gal4
	@echo "\n\e[35m\e[1m== Node Omega restart =======================================================\e[0m"
	bin/nodeomega _test/nodeomega_gal4
	@echo "\n\e[35m\e[1m== Node MutSel run ==========================================================\e[0m"
	bin/nodemutsel --ncat 3 -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Node MutSel restart ======================================================\e[0m"
	bin/nodemutsel _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Node MutSel read =========================================================\e[0m"
	bin/readnodemutsel --ss _test/nodemutsel_gal4

# @make --no-print-directory run-multigeneglobom-test
.PHONY: run-multigeneglobom-test
run-multigeneglobom-test: all
	@echo "\n\e[35m\e[1m== Multigene Single Omega ===================================================\e[0m"
	cd data/small_multigene && mpirun -np 2 ../../bin/multigeneglobom -t tree.nwk -a verysmall.list  -u ${POINTS} tmp

.PHONY: testpr
test:
	@make --no-print-directory run-unit-tests
	@make --no-print-directory run-unit-tests-mpi
	@make --no-print-directory run-app-tests

.PHONY: aamutsel
aamutsel: bin
	@cd bin ; make --no-print-directory -j8 aamutsel readaamutsel
	@rm -rf _aamutsel
	@mkdir _aamutsel
	bin/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 _aamutsel/gal4
	bin/aamutsel _aamutsel/gal4
	bin/readaamutsel _aamutsel/gal4
	bin/readaamutsel --ss _aamutsel/gal4
	bin/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 -p _aamutsel/gal4_poly
	bin/aamutsel _aamutsel/gal4_poly
	bin/readaamutsel _aamutsel/gal4_poly
	bin/readaamutsel --ss _aamutsel/gal4_poly

.PHONY: mutselomega
mutselomega: tiny
	@cd bin ; make --no-print-directory -j8 mutselomega readmutselomega
	@rm -rf _mutselomega
	@mkdir _mutselomega
	bin/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --ncat 30 -u 30 _mutselomega/mutsel_bglobin
	bin/mutselomega _mutselomega/mutsel_bglobin
	bin/readmutselomega -b 10 --ss _mutselomega/mutsel_bglobin
	bin/readmutselomega -b 10 --omega_0 _mutselomega/mutsel_bglobin
	bin/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --freeomega --omegancat 3 -u 30 --flatfitness _mutselomega/MGM3_bglobin
	bin/mutselomega _mutselomega/MGM3_bglobin
	bin/readmutselomega -b 10 _mutselomega/MGM3_bglobin
	bin/readmutselomega -b 10 _mutselomega/mutsel_bglobin --chain_omega _mutselomega/MGM3_bglobin
	bin/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --profiles _mutselomega/mutsel_bglobin.siteprofiles --omegaarray data/bglobin/omegaarray.csv -u 30 _mutselomega/clamped_bglobin
	bin/mutselomega _mutselomega/clamped_bglobin
	bin/readmutselomega -b 10 _mutselomega/clamped_bglobin
	bin/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --freeomega --omegancat 3 --ncat 30 -u 30 _mutselomega/mutselM3_bglobin
	bin/mutselomega _mutselomega/mutselM3_bglobin
	bin/readmutselomega -b 10 _mutselomega/mutselM3_bglobin

.PHONY: DM5
DM5: bin
	@cd bin ; make --no-print-directory -j8 mutseldm5 readmutseldm5
	@rm -rf _mutseldm5
	@mkdir _mutseldm5
	bin/mutseldm5 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre  --freeomega --omegancat 10 --flatfitness --fixp0 --p0 0.0 -u 30 --hypermean_threshold 0.0 --hyperinvshape_threshold 10.0 _mutseldm5/MGM3_bglobin
	bin/readmutseldm5 _mutseldm5/MGM3_bglobin
	bin/mutseldm5 _mutseldm5/MGM3_bglobin
	bin/readmutseldm5 --ss _mutseldm5/MGM3_bglobin
	bin/mutseldm5 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre  --omegashift 1.0 --freeomega --omegancat 10 --ncat 30 -u 30 _mutseldm5/bglobin
	bin/readmutseldm5 _mutseldm5/bglobin
	bin/mutseldm5 _mutseldm5/bglobin

.PHONY: dated
dated: tiny
	@cd bin ; make --no-print-directory -j8 nodemutsel
	@rm -rf _dated
	@mkdir _dated
	bin/nodemutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --ncat 3 -u ${POINTS} _dated/node_gal4
	bin/nodemutsel _dated/node_gal4
	bin/nodemutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --ncat 3 -u ${POINTS} -p _dated/node_poly_gal4
	bin/nodemutsel _dated/node_poly_gal4

.PHONY: traits
traits: tiny
	@cd bin ; make --no-print-directory -j8 nodetraits readnodetraits
	@rm -rf _traits
	@mkdir _traits
	python3 utils/neutrality_index.py --tree data/body_size/mammals.male.tree --traits data/body_size/mammals.male.traits.tsv --var_within data/body_size/mammals.male.var_within.tsv --output _traits/mammals.male.ML.tsv
	bin/nodetraits -t data/body_size/mammals.male.tree --traitsfile data/body_size/mammals.male.traits.tsv -u 100 _traits/mammals.male
	bin/nodetraits _traits/mammals.male
	bin/readnodetraits -b 50 --newick _traits/mammals.male
	bin/readnodetraits -b 50 --cov _traits/mammals.male
	python3 utils/ratio_pvalue.py --burn_in 50 --inference _traits/mammals.male.trace --var_within data/body_size/mammals.male.var_within.tsv --output _traits/mammals.male.Bayes.tsv

.PHONY: diffseldsparse
diffseldsparse: all
		@rm -f delme*.*
		bin/diffseldsparse -a data/besnard/cyp_small.phy -t data/besnard/cyp_coding.Chrysithr_root.nhx -e 1 -u 3 tmp
