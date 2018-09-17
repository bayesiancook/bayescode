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
	@clang-format -i sources/*.cpp sources/*.hpp