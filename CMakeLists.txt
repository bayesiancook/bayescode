cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(bayescode CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter") # normal mode
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src")
include_directories("src/components")
include_directories("src/tree")
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
    src/lib/PhyloProcess.cpp
    src/lib/Random.cpp
    src/lib/RandomGeneSample.cpp
    src/lib/Sample.cpp
    src/lib/SequenceAlignment.cpp
    src/lib/StateSpace.cpp
    src/lib/SubMatrix.cpp
    src/lib/TaxonSet.cpp
    src/tree/tree-implem.cpp
  )
add_library (bayescode_lib STATIC ${BAYESCODE_LIB})
target_compile_features(bayescode_lib PUBLIC ${CXX_FEATURES})

# Globom
add_executable(globom "src/SingleOmega.cpp")
target_link_libraries(globom bayescode_lib)
target_compile_features(globom PRIVATE ${CXX_FEATURES})

# Diffsel double sparse
add_executable(diffseldsparse "src/DiffSelDoublySparse.cpp")
target_link_libraries(diffseldsparse bayescode_lib)
target_compile_features(diffseldsparse PRIVATE ${CXX_FEATURES})

#TODO readsingleomega and readdiffseldoublesparse