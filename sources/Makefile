#To compile with gsl, reactivate the following line (remove '#')
#USE_GSL=-DGSL

#To statically build the executables, uncomment the following line
#STATIC=-static

ifdef USE_GSL
#GSLLIB=-lgsl -lgslcblas
else
GSLLIB=
endif

LIBS=$(GSLLIB)
#LIBS=$(GSLLIB) -lboost_regex

CC=mpic++
SYSLIB=
INCLUDES=
#CPPFLAGS= $(USE_GSL) $(INCLUDES)
CPPFLAGS= -std=c++11 -g -Wall -O3 -Wno-int-in-bool-context $(USE_GSL) $(INCLUDES)
LDFLAGS=$(STATIC)
INSTALL_DIR=
INSTALL_LIB=
SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp AAMutSelOmegaCodonSubMatrix.cpp GTRSubMatrix.cpp AASubSelSubMatrix.cpp AAMutSelSubMatrix.cpp T92SubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp linalg.cpp cdf.cpp Chain.cpp MultiGeneChain.cpp Sample.cpp MultiGeneSample.cpp MPIBuffer.cpp MultiGeneMPIModule.cpp CodonM2aModel.cpp MultiGeneCodonM2aModel.cpp AAMutSelM2aModel.cpp MultiGeneAAMutSelM2aModel.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= ppredtest multigenesiteom siteom readmultigenediffseldsparse readmultigeneaamutselddp multigenecodonm2a sparseaamutsel sparseaasubsel sparsecatgtr simucodonm2a readcodonm2a aamutselsparse diffseldsparse readdiffseldsparse readglobom multigenediffseldsparse multigenediffselsparse multigenecondom multigeneaamutselddp multigenediffselsparse codonm2a aamutselddp readaamutselddp diffselsparse readdiffselsparse diffsel multigenediffsel readdiffsel multigenebranchom readmultigenebranchom multigeneglobom globom readglobom readmultigeneglobom readmultigenecodonm2a
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

# If we are on a windows platform, executables are .exe files
# [S.Bigaret] we assume here that we are compiling w/ MinGW, not Cygwin
EXEEXT=
uname := $(shell uname -s 2>/dev/null)
ifeq (MINGW,$(findstring MINGW,$(uname)))
	EXEEXT=.exe
	PROGS=$(addsuffix $(EXEEXT), $(ALL))
endif

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

subsample$(EXEEXT): $(PROGSDIR)/subsample$(EXEEXT)
$(PROGSDIR)/subsample$(EXEEXT): RandomGeneSample.o $(OBJS)
	$(CC) RandomGeneSample.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

chi2$(EXEEXT): $(PROGSDIR)/chi2$(EXEEXT)
$(PROGSDIR)/chi2$(EXEEXT): Chi2.o $(OBJS)
	$(CC) Chi2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ppredtest$(EXEEXT): $(PROGSDIR)/ppredtest$(EXEEXT)
$(PROGSDIR)/ppredtest$(EXEEXT): PostPredTest.o $(OBJS)
	$(CC) PostPredTest.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

sparseaamutsel$(EXEEXT): $(PROGSDIR)/sparseaamutsel$(EXEEXT)
$(PROGSDIR)/sparseaamutsel$(EXEEXT): SparseAAMutSel.o $(OBJS)
	$(CC) SparseAAMutSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

sparseaasubsel$(EXEEXT): $(PROGSDIR)/sparseaasubsel$(EXEEXT)
$(PROGSDIR)/sparseaasubsel$(EXEEXT): SparseAASubSel.o $(OBJS)
	$(CC) SparseAASubSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

sparsecatgtr$(EXEEXT): $(PROGSDIR)/sparsecatgtr$(EXEEXT)
$(PROGSDIR)/sparsecatgtr$(EXEEXT): SparseCATGTR.o $(OBJS)
	$(CC) SparseCATGTR.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselsparse$(EXEEXT): $(PROGSDIR)/aamutselsparse$(EXEEXT)
$(PROGSDIR)/aamutselsparse$(EXEEXT): AAMutSelSparseOmega.o $(OBJS)
	$(CC) AAMutSelSparseOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readaamutselddp$(EXEEXT): $(PROGSDIR)/readaamutselddp$(EXEEXT)
$(PROGSDIR)/readaamutselddp$(EXEEXT): ReadAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) ReadAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselddp$(EXEEXT): $(PROGSDIR)/aamutselddp$(EXEEXT)
$(PROGSDIR)/aamutselddp$(EXEEXT): AAMutSelDSBDPOmega.o $(OBJS)
	$(CC) AAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselddp$(EXEEXT): $(PROGSDIR)/multigeneaamutselddp$(EXEEXT)
$(PROGSDIR)/multigeneaamutselddp$(EXEEXT): MultiGeneAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) MultiGeneAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselm2a$(EXEEXT): $(PROGSDIR)/aamuselm2a$(EXEEXT)
$(PROGSDIR)/aamutselm2a$(EXEEXT): AAMutSelM2a.o $(OBJS)
	$(CC) AAMutSelM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simucodonm2a$(EXEEXT): $(PROGSDIR)/simucodonm2a$(EXEEXT)
$(PROGSDIR)/simucodonm2a$(EXEEXT): SimuCodonM2a.o $(OBJS)
	$(CC) SimuCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm2a$(EXEEXT): $(PROGSDIR)/codonm2a$(EXEEXT)
$(PROGSDIR)/codonm2a$(EXEEXT): CodonM2a.o $(OBJS)
	$(CC) CodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm2a$(EXEEXT): $(PROGSDIR)/multigenecodonm2a$(EXEEXT)
$(PROGSDIR)/multigenecodonm2a$(EXEEXT): MultiGeneCodonM2a.o $(OBJS)
	$(CC) MultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselm2a$(EXEEXT): $(PROGSDIR)/multigeneaamutselm2a$(EXEEXT)
$(PROGSDIR)/multigeneaamutselm2a$(EXEEXT): MultiGeneAAMutSelM2a.o $(OBJS)
	$(CC) MultiGeneAAMutSelM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneglobom$(EXEEXT): $(PROGSDIR)/multigeneglobom$(EXEEXT)
$(PROGSDIR)/multigeneglobom$(EXEEXT): MultiGeneSingleOmega.o $(OBJS)
	$(CC) MultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

dgamom$(EXEEXT): $(PROGSDIR)/dgamom$(EXEEXT)
$(PROGSDIR)/dgamom$(EXEEXT): DiscGammaSiteOmega.o $(OBJS)
	$(CC) DiscGammaSiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

dpom$(EXEEXT): $(PROGSDIR)/dpom$(EXEEXT)
$(PROGSDIR)/dpom$(EXEEXT): DPOmega.o $(OBJS)
	$(CC) DPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenediffseldsparse$(EXEEXT): $(PROGSDIR)/readmultigenediffseldsparse$(EXEEXT)
$(PROGSDIR)/readmultigenediffseldsparse$(EXEEXT): ReadMultiGeneDiffSelDoublySparse.o $(OBJS)
	$(CC) ReadMultiGeneDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigeneaamutselddp$(EXEEXT): $(PROGSDIR)/readmultigeneaamutselddp$(EXEEXT)
$(PROGSDIR)/readmultigeneaamutselddp$(EXEEXT): ReadMultiGeneAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) ReadMultiGeneAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigeneglobom$(EXEEXT): $(PROGSDIR)/readmultigeneglobom$(EXEEXT)
$(PROGSDIR)/readmultigeneglobom$(EXEEXT): ReadMultiGeneSingleOmega.o $(OBJS)
	$(CC) ReadMultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenecodonm2a$(EXEEXT): $(PROGSDIR)/readmultigenecodonm2a$(EXEEXT)
$(PROGSDIR)/readmultigenecodonm2a$(EXEEXT): ReadMultiGeneCodonM2a.o $(OBJS)
	$(CC) ReadMultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenebranchom$(EXEEXT): $(PROGSDIR)/readmultigenebranchom$(EXEEXT)
$(PROGSDIR)/readmultigenebranchom$(EXEEXT): ReadMultiGeneBranchOmega.o $(OBJS)
	$(CC) ReadMultiGeneBranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenebranchom$(EXEEXT): $(PROGSDIR)/multigenebranchom$(EXEEXT)
$(PROGSDIR)/multigenebranchom$(EXEEXT): MultiGeneBranchOmega.o $(OBJS)
	$(CC) MultiGeneBranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecondom$(EXEEXT): $(PROGSDIR)/multigenecondom$(EXEEXT)
$(PROGSDIR)/multigenecondom$(EXEEXT): MultiGeneConditionOmega.o $(OBJS)
	$(CC) MultiGeneConditionOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

globom$(EXEEXT): $(PROGSDIR)/globom$(EXEEXT)
$(PROGSDIR)/globom$(EXEEXT): SingleOmega.o $(OBJS)
	$(CC) SingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcodonm2a$(EXEEXT): $(PROGSDIR)/readcodonm2a$(EXEEXT)
$(PROGSDIR)/readcodonm2a$(EXEEXT): ReadCodonM2a.o $(OBJS)
	$(CC) ReadCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readglobom$(EXEEXT): $(PROGSDIR)/readglobom$(EXEEXT)
$(PROGSDIR)/readglobom$(EXEEXT): ReadSingleOmega.o $(OBJS)
	$(CC) ReadSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenesiteom$(EXEEXT): $(PROGSDIR)/multigenesiteom$(EXEEXT)
$(PROGSDIR)/multigenesiteom$(EXEEXT): MultiGeneSiteOmega.o $(OBJS)
	$(CC) MultiGeneSiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

siteom$(EXEEXT): $(PROGSDIR)/siteom$(EXEEXT)
$(PROGSDIR)/siteom$(EXEEXT): SiteOmega.o $(OBJS)
	$(CC) SiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiffseldsparse$(EXEEXT): $(PROGSDIR)/readdiffseldsparse$(EXEEXT)
$(PROGSDIR)/readdiffseldsparse$(EXEEXT): ReadDiffSelDoublySparse.o $(OBJS)
	$(CC) ReadDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiffselsparse$(EXEEXT): $(PROGSDIR)/readdiffselsparse$(EXEEXT)
$(PROGSDIR)/readdiffselsparse$(EXEEXT): ReadDiffSelSparse.o $(OBJS)
	$(CC) ReadDiffSelSparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiffsel$(EXEEXT): $(PROGSDIR)/readdiffsel$(EXEEXT)
$(PROGSDIR)/readdiffsel$(EXEEXT): ReadDiffSel.o $(OBJS)
	$(CC) ReadDiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffseldsparse$(EXEEXT): $(PROGSDIR)/diffseldsparse$(EXEEXT)
$(PROGSDIR)/diffseldsparse$(EXEEXT): DiffSelDoublySparse.o $(OBJS)
	$(CC) DiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffselsparse$(EXEEXT): $(PROGSDIR)/diffselsparse$(EXEEXT)
$(PROGSDIR)/diffselsparse$(EXEEXT): DiffSelSparse.o $(OBJS)
	$(CC) DiffSelSparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediffselsparse$(EXEEXT): $(PROGSDIR)/multigenediffselsparse$(EXEEXT)
$(PROGSDIR)/multigenediffselsparse$(EXEEXT): MultiGeneDiffSelSparse.o $(OBJS)
	$(CC) MultiGeneDiffSelSparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediffsel$(EXEEXT): $(PROGSDIR)/multigenediffsel$(EXEEXT)
$(PROGSDIR)/multigenediffsel$(EXEEXT): MultiGeneDiffSel.o $(OBJS)
	$(CC) MultiGeneDiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediffseldsparse$(EXEEXT): $(PROGSDIR)/multigenediffseldsparse$(EXEEXT)
$(PROGSDIR)/multigenediffseldsparse$(EXEEXT): MultiGeneDiffSelDoublySparse.o $(OBJS)
	$(CC) MultiGeneDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffsel$(EXEEXT): $(PROGSDIR)/diffsel$(EXEEXT)
$(PROGSDIR)/diffsel$(EXEEXT): DiffSel.o $(OBJS)
	$(CC) DiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

.PHONY: format
format:
	clang-format -i *.hpp *.cpp
