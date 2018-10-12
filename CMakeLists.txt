cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(bayescode CXX)

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter -Wfatal-errors -Wno-int-in-bool-context") # normal mode

# Sanitize: runtime checks for various forms of undefined or suspicious behavior (clang++).
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -g")
# For [-fsanitize] use one of the following:
# -fsanitize=address: AddressSanitizer, a memory error detector.
# -fsanitize=thread: ThreadSanitizer, a data race detector.
# -fsanitize=memory: MemorySanitizer, a detector of uninitialized reads. Requires instrumentation of all program code.
# -fsanitize=undefined: UndefinedBehaviorSanitizer, a fast and compatible undefined behavior checker.
# -fsanitize=dataflow: DataFlowSanitizer, a general data flow analysis.
# -fsanitize=cfi: control flow integrity checks. Requires -flto.
# -fsanitize=safe-stack: safe stack protection against stack-based memory corruption errors.
# Use [-g] to have file names and line numbers included in warning messages.

# Read AAMutSel, but site profiles is not working yet
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src")
include_directories("src/lib")
include_directories("utils")

# Static lib
set(BAYESCODE_LIB
    src/lib/AAMutSelOmegaCodonSubMatrix.cpp
    src/lib/BranchSitePath.cpp
    src/lib/Chain.cpp
    src/lib/Chrono.cpp
    src/lib/CodonSequenceAlignment.cpp
    src/lib/CodonStateSpace.cpp
    src/lib/CodonSubMatrix.cpp
    src/lib/GTRSubMatrix.cpp
    src/lib/linalg.cpp
    src/lib/MPIBuffer.cpp
    src/lib/MultiGeneChain.cpp
    src/lib/MultiGeneMPIModule.cpp
    src/lib/PhyloProcess.cpp
    src/lib/Random.cpp
    src/lib/RandomGeneSample.cpp
    src/lib/Sample.cpp
    src/lib/SequenceAlignment.cpp
    src/lib/StateSpace.cpp
    src/lib/SubMatrix.cpp
    src/lib/TaxonSet.cpp
    src/CodonM2aModel.cpp
    src/tree/implem.cpp
  )
add_library (bayescode_lib STATIC ${BAYESCODE_LIB})
target_compile_features(bayescode_lib PUBLIC ${CXX_FEATURES})

# Single omega (globom)
add_executable(globom "src/SingleOmega.cpp")
target_link_libraries(globom bayescode_lib)
target_compile_features(globom PRIVATE ${CXX_FEATURES})

# Codon M2a
add_executable(codonm2a "src/CodonM2a.cpp")
target_link_libraries(codonm2a bayescode_lib)
target_compile_features(codonm2a PRIVATE ${CXX_FEATURES})

# Read Codon M2a
add_executable(readcodonm2a "src/ReadCodonM2a.cpp")
target_link_libraries(readcodonm2a bayescode_lib)
target_compile_features(readcodonm2a PRIVATE ${CXX_FEATURES})

# Mutation-Selection (aamutsel)
add_executable(aamutsel "src/AAMutSelDSBDPOmega.cpp")
target_link_libraries(aamutsel bayescode_lib)
target_compile_features(aamutsel PRIVATE ${CXX_FEATURES})

# Read Mutation-Selection (aamutsel)
add_executable(readaamutsel "src/ReadAAMutSelDSBDPOmega.cpp")
target_link_libraries(readaamutsel bayescode_lib)
target_compile_features(readaamutsel PRIVATE ${CXX_FEATURES})

# Diffsel double sparse (diffseldsparse)
add_executable(diffseldsparse "src/DiffSelDoublySparse.cpp")
target_link_libraries(diffseldsparse bayescode_lib)
target_compile_features(diffseldsparse PRIVATE ${CXX_FEATURES})

# Read diffsel double sparse (diffseldsparse)
add_executable(readdiffseldsparse "src/ReadDiffSelDoublySparse.cpp")
target_link_libraries(readdiffseldsparse bayescode_lib)
target_compile_features(readdiffseldsparse PRIVATE ${CXX_FEATURES})

# Multi-gene single omega (multigeneglobom)
add_executable(multigeneglobom "src/MultiGeneSingleOmega.cpp")
target_link_libraries(multigeneglobom bayescode_lib ${MPI_LIBRARIES})
target_compile_features(multigeneglobom PRIVATE ${CXX_FEATURES})

# New Multi-gene single omega (multigeneglobomnew)
add_executable(multigeneglobomnew "src/MultiGeneSingleOmegaNew.cpp")
target_link_libraries(multigeneglobomnew bayescode_lib ${MPI_LIBRARIES})
target_compile_features(multigeneglobomnew PRIVATE ${CXX_FEATURES})

#TODO readdiffseldoublesparse
