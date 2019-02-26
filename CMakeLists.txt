cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(bayescode CXX)

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Compilation options
option(COVERAGE_MODE "For coverage mode using g++ " OFF) #OFF by default
if(COVERAGE_MODE)
    set(CMAKE_CXX_FLAGS "-O0 --std=c++14 -fprofile-arcs -ftest-coverage") # coverage mode
else()
    set(CMAKE_CXX_FLAGS "-O3 --std=c++14 -Wall -Wextra -Wno-unused-parameter -Wfatal-errors") # normal mode
endif(COVERAGE_MODE)
unset(COVERAGE_MODE CACHE)

include_directories("src")
include_directories("src/lib")
include_directories("utils")

# Static libs
add_library (tree_lib STATIC "src/tree/implem.cpp")
add_library (mpicomponent_lib STATIC "src/mpi_components/Process.cpp")

set(BAYESCODE_LIB
    src/lib/AAMutSelOmegaCodonSubMatrix.cpp
    src/lib/BranchSitePath.cpp
    src/lib/Chain.cpp
    src/lib/Chrono.cpp
    src/lib/CodonSequenceAlignment.cpp
    src/lib/CodonStateSpace.cpp
    src/lib/CodonSubMatrix.cpp
    src/lib/GTRSubMatrix.cpp
    src/lib/PhyloProcess.cpp
    src/lib/Random.cpp
    src/lib/RandomGeneSample.cpp
    src/lib/Sample.cpp
    src/lib/SequenceAlignment.cpp
    src/lib/StateSpace.cpp
    src/lib/SubMatrix.cpp
    src/lib/TaxonSet.cpp
    src/lib/PolyData.cpp
    src/lib/PolyProcess.cpp
    src/lib/PoissonRandomField.cpp
  )
add_library (bayescode_lib STATIC ${BAYESCODE_LIB})

# Single omega (globom)
add_executable(globom "src/SingleOmega.cpp")
target_link_libraries(globom bayescode_lib tree_lib)

# Read single omega (readglobom)
add_executable(readglobom "src/ReadSingleOmega.cpp")
target_link_libraries(readglobom bayescode_lib tree_lib)

# Codon M2a
add_executable(codonm2a "src/CodonM2a.cpp" "src/CodonM2aModel.cpp")
target_link_libraries(codonm2a bayescode_lib tree_lib)

# Read Codon M2a
add_executable(readcodonm2a "src/ReadCodonM2a.cpp" "src/CodonM2aModel.cpp")
target_link_libraries(readcodonm2a bayescode_lib tree_lib)

# Mutation-Selection (aamutsel)
add_executable(aamutsel "src/AAMutSelDSBDPOmega.cpp")
target_link_libraries(aamutsel bayescode_lib tree_lib)

# Read Mutation-Selection (readaamutsel)
add_executable(readaamutsel "src/ReadAAMutSelDSBDPOmega.cpp")
target_link_libraries(readaamutsel bayescode_lib tree_lib)

# Mutation-Selection (mutselomega)
add_executable(mutselomega "src/AAMutSelMultipleOmega.cpp")
target_link_libraries(mutselomega bayescode_lib tree_lib)

# Read Mutation-Selection (readmutselomega)
add_executable(readmutselomega "src/ReadAAMutSelMultipleOmega.cpp")
target_link_libraries(readmutselomega bayescode_lib tree_lib)

# Diffsel double sparse (diffseldsparse)
add_executable(diffseldsparse "src/DiffSelDoublySparse.cpp")
target_link_libraries(diffseldsparse bayescode_lib tree_lib)

# Read diffsel double sparse (readdiffseldsparse)
# add_executable(readdiffseldsparse "src/ReadDiffSelDoublySparse.cpp")
# target_link_libraries(readdiffseldsparse bayescode_lib tree_lib)

# Multi-gene single omega (multigeneglobom)
add_executable(multigeneglobom "src/MultiGeneSingleOmega.cpp")
target_link_libraries(multigeneglobom bayescode_lib tree_lib ${MPI_LIBRARIES} mpicomponent_lib)

# Multi Gene Codon M2a
add_executable(multigenecodonm2a "src/MultiGeneCodonM2a.cpp" "src/CodonM2aModel.cpp")
target_link_libraries(multigenecodonm2a bayescode_lib tree_lib ${MPI_LIBRARIES} mpicomponent_lib)

# tests
add_executable(all_tests "src/all_tests.cpp")
target_link_libraries(all_tests ${MPI_LIBRARIES} mpicomponent_lib)

add_executable(tree_test "src/tree/test.cpp")
target_link_libraries(tree_test tree_lib)

add_executable(mpi_par_test "src/mpi_components/test/mpi_test.cpp")
target_link_libraries(mpi_par_test ${MPI_LIBRARIES} mpicomponent_lib)
