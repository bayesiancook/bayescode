project('bayescode', ['cpp','c'],
	version: '0.0.1',
	default_options : [
	  'buildtype=release',
	  'cpp_std=c++11'
	],
	meson_version: '>= 0.49',
	license: 'unknown')

#To compile with gsl, reactivate the following line (remove '#')

cpp = meson.get_compiler('cpp')
cc = meson.get_compiler('c')

warn_flags = ['-Wall'] # ['-Wextra','-Wno-sign-compare','-Wno-maybe-uninitialized','-Woverloaded-virtual','-Wstrict-aliasing','-Wno-unknown-pragmas', '-fdiagnostics-show-template-tree']

add_project_arguments(cpp.get_supported_arguments(warn_flags), language : 'cpp')

if get_option('buildtype') == 'release' or get_option('buildtype') == 'plain'
  opt_flags = ['-O3','-funroll-loops','-fno-math-errno','-fno-signed-zeros']
  add_project_arguments(['-DNDEBUG']+cpp.get_supported_arguments(opt_flags), language : 'cpp')
endif

mpi = dependency('mpi', language: 'cpp')

srcs= [ 'BranchSitePath.cpp', 'Chrono.cpp', 'CodonSequenceAlignment.cpp', 'CodonStateSpace.cpp', 'CodonSubMatrix.cpp', 'AAMutSelOmegaCodonSubMatrix.cpp', 'GTRSubMatrix.cpp', 'AASubSelSubMatrix.cpp', 'AAMutSelSubMatrix.cpp', 'T92SubMatrix.cpp', 'PhyloProcess.cpp', 'Random.cpp', 'SequenceAlignment.cpp', 'StateSpace.cpp', 'SubMatrix.cpp', 'TaxonSet.cpp', 'Tree.cpp', 'linalg.cpp', 'cdf.cpp', 'Chain.cpp', 'MultiGeneChain.cpp', 'Sample.cpp', 'MultiGeneSample.cpp', 'MPIBuffer.cpp', 'MultiGeneMPIModule.cpp', 'CodonM2aModel.cpp', 'MultiGeneCodonM2aModel.cpp']

libcommon = static_library('common',
                           srcs,
                           include_directories: [],
                           dependencies: [mpi])

progs =   { 'globom': ['SingleOmega.cpp'],
            'readglobom': ['ReadSingleOmega.cpp'],
            'multigeneglobom': ['MultiGeneSingleOmega.cpp'],
            'readmultigeneglobom': ['ReadMultiGeneSingleOmega.cpp'],
            'codonm2a': ['CodonM2a.cpp'],
            'readcodonm2a': ['ReadCodonM2a.cpp'],
            'simucodonm2a': ['SimuCodonM2a.cpp'],
            'multigenecodonm2a': ['MultiGeneCodonM2a.cpp'],
            'readmultigenecodonm2a': ['ReadMultiGeneCodonM2a.cpp'],
            'fastreadmultigenecodonm2a': ['FastReadMultiGeneCodonM2a.cpp'],
            'aamutselddp': ['AAMutSelDSBDPOmega.cpp'],
            'readaamutselddp': ['ReadAAMutSelDSBDPOmega.cpp'],
            'multigeneaamutselddp': ['MultiGeneAAMutSelDSBDPOmega.cpp'],
            'readmultigeneaamutselddp': ['ReadMultiGeneAAMutSelDSBDPOmega.cpp'],
            'diffsel': ['DiffSel.cpp'],
            'readdiffsel': ['ReadDiffSel.cpp'],
            'multigenediffsel': ['MultiGeneDiffSel.cpp'],
            'diffseldsparse': ['DiffSelDoublySparse.cpp'],
            'readdiffseldsparse': ['ReadDiffSelDoublySparse.cpp'],
            'multigenediffseldsparse': ['MultiGeneDiffSelDoublySparse.cpp'],
            'readmultigenediffseldsparse': ['ReadMultiGeneDiffSelDoublySparse.cpp'],
            'multigenebranchom': ['MultiGeneConditionOmega.cpp'],
            'readmultigenebranchom': ['ReadMultiGeneConditionOmega.cpp'],
            'multigenesparsebranchom': ['MultiGeneSparseConditionOmega.cpp'],
            'readmultigenesparsebranchom': ['ReadMultiGeneSparseConditionOmega.cpp'],
            'ppredtest': ['PostPredTest.cpp'],
            'multigenesiteom': ['MultiGeneSiteOmega.cpp'],
            'siteom': ['SiteOmega.cpp']
          }

foreach name,sources: progs
  exe = executable(name,
                   sources,
                   dependencies: [mpi],
                   link_with: [libcommon],
                   install: true)
endforeach


