cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(bayescode CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter -Wno-inconsistent-missing-override") # normal mode
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("sources")

# static lib
set(BAYESCODE_LIB
    sources/AAMutSelOmegaCodonSubMatrix.cpp
    sources/BranchSitePath.cpp
    sources/Chain.cpp
    sources/Chrono.cpp
    sources/CodonSequenceAlignment.cpp
    sources/CodonStateSpace.cpp
    sources/CodonSubMatrix.cpp
    sources/GTRSubMatrix.cpp
    sources/linalg.cpp
    sources/MPIBuffer.cpp
    sources/PhyloProcess.cpp
    sources/Random.cpp
    sources/RandomGeneSample.cpp
    sources/Sample.cpp
    sources/SequenceAlignment.cpp
    sources/StateSpace.cpp
    sources/SubMatrix.cpp
    sources/TaxonSet.cpp
    sources/tree-implem.cpp
  )
add_library (bayescode_lib STATIC ${BAYESCODE_LIB})
target_compile_features(bayescode_lib PUBLIC ${CXX_FEATURES})

# Globom
add_executable(globom "sources/SingleOmega.cpp")
target_link_libraries(globom bayescode_lib)
target_compile_features(globom PRIVATE ${CXX_FEATURES})

# Globom
add_executable(diffseldsparse "sources/DiffSelDoublySparse.cpp")
target_link_libraries(diffseldsparse bayescode_lib)
target_compile_features(diffseldsparse PRIVATE ${CXX_FEATURES})

#TODO readsingleomega and readdiffseldoublesparse