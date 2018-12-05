all:
	@cd sources && make -j8 --no-print-directory

.PHONY: format
format:
	clang-format -i sources/*.hpp sources/*.cpp

.PHONY: test-mgdiffsel
test-mgdiffsel: all
	cd data/small_multigene/ && mpirun -np 3 ../multigenediffseldsparse -t tree.nwk.annotated -d verysmall.list -ncond 2 -x 1 10 tmp