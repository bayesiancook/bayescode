# COMPILATION
# Requires: cmake 3.1.0 or better
.PHONY:all
all: cmake
	@cd _build ; make --no-print-directory -j8

.PHONY: cmake
cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

.PHONY:clean
clean:
	@rm -rf _build

# CODE QUALITY
# Requires: clang-format
.PHONY: format
format:
	@clang-format -i `find -name *.*pp`

# Run AaMutSel
.PHONY: aamutsel
aamutsel: all
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 data/polymorphism/gal4_1;
	_build/readaamutsel --om -b 0 -e 1 -u 10 data/polymorphism/gal4_1;
	_build/readaamutsel --ss -b 0 -e 1 -u 10 data/polymorphism/gal4_1;