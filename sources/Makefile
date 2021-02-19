CC=mpic++
CPPFLAGS= -std=c++11 -Wall -O3
LDFLAGS= -O3
SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp AAMutSelOmegaCodonSubMatrix.cpp GTRSubMatrix.cpp AASubSelSubMatrix.cpp AAMutSelSubMatrix.cpp T92SubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp linalg.cpp cdf.cpp Chain.cpp MultiGeneChain.cpp Sample.cpp MultiGeneSample.cpp MPIBuffer.cpp MultiGeneMPIModule.cpp CovMatrix.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= epsrobustbranchom readepsrobustbranchom branchom readbranchom coevol readcoevol brownclock globom codonm2a readaamutselsparse readaamutselddp aamutselsparsem9 multigeneaamutselsparsem9 codonm9 multigenecodonm9 codonm8 multigeneaamutselsparse aamutselsparse aamutselgdp aamutseldp fastempselac selempfreq fastselac multigeneaamutselm2 aamutselm2 discselac multigenediscselac empselac selac multigeneselac globom readglobom multigeneglobom readmultigeneglobom codonm2a readcodonm2a simucodonm2a multigenecodonm2a readmultigenecodonm2a fastreadmultigenecodonm2a aamutselddp readaamutselddp multigeneaamutselddp readmultigeneaamutselddp diffsel readdiffsel multigenediffsel diffseldsparse readdiffseldsparse multigenediffseldsparse readmultigenediffseldsparse multigenebranchom readmultigenebranchom multigenesparsebranchom readmultigenesparsebranchom ppredtest multigenesiteom siteom 
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

epsrobustbranchom$(EXEEXT): $(PROGSDIR)/epsrobustbranchom$(EXEEXT)
$(PROGSDIR)/epsrobustbranchom$(EXEEXT): EpsRobustBranchOmega.o $(OBJS)
	$(CC) EpsRobustBranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readepsrobustbranchom$(EXEEXT): $(PROGSDIR)/readepsrobustbranchom$(EXEEXT)
$(PROGSDIR)/readepsrobustbranchom$(EXEEXT): ReadEpsRobustBranchOmega.o $(OBJS)
	$(CC) ReadEpsRobustBranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

branchom$(EXEEXT): $(PROGSDIR)/branchom$(EXEEXT)
$(PROGSDIR)/branchom$(EXEEXT): BranchOmega.o $(OBJS)
	$(CC) BranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readbranchom$(EXEEXT): $(PROGSDIR)/readbranchom$(EXEEXT)
$(PROGSDIR)/readbranchom$(EXEEXT): ReadBranchOmega.o $(OBJS)
	$(CC) ReadBranchOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

coevol$(EXEEXT): $(PROGSDIR)/coevol$(EXEEXT)
$(PROGSDIR)/coevol$(EXEEXT): Coevol.o $(OBJS)
	$(CC) Coevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcoevol$(EXEEXT): $(PROGSDIR)/readcoevol$(EXEEXT)
$(PROGSDIR)/readcoevol$(EXEEXT): ReadCoevol.o $(OBJS)
	$(CC) ReadCoevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

brownclock$(EXEEXT): $(PROGSDIR)/brownclock$(EXEEXT)
$(PROGSDIR)/brownclock$(EXEEXT): BrownianClock.o $(OBJS)
	$(CC) BrownianClock.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

globom$(EXEEXT): $(PROGSDIR)/globom$(EXEEXT)
$(PROGSDIR)/globom$(EXEEXT): SingleOmega.o $(OBJS)
	$(CC) SingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readglobom$(EXEEXT): $(PROGSDIR)/readglobom$(EXEEXT)
$(PROGSDIR)/readglobom$(EXEEXT): ReadSingleOmega.o $(OBJS)
	$(CC) ReadSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneglobom$(EXEEXT): $(PROGSDIR)/multigeneglobom$(EXEEXT)
$(PROGSDIR)/multigeneglobom$(EXEEXT): MultiGeneSingleOmega.o $(OBJS)
	$(CC) MultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigeneglobom$(EXEEXT): $(PROGSDIR)/readmultigeneglobom$(EXEEXT)
$(PROGSDIR)/readmultigeneglobom$(EXEEXT): ReadMultiGeneSingleOmega.o $(OBJS)
	$(CC) ReadMultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm9$(EXEEXT): $(PROGSDIR)/codonm9$(EXEEXT)
$(PROGSDIR)/codonm9$(EXEEXT): CodonM9.o $(OBJS)
	$(CC) CodonM9.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm8$(EXEEXT): $(PROGSDIR)/codonm8$(EXEEXT)
$(PROGSDIR)/codonm8$(EXEEXT): CodonM8.o $(OBJS)
	$(CC) CodonM8.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm9$(EXEEXT): $(PROGSDIR)/multigenecodonm9$(EXEEXT)
$(PROGSDIR)/multigenecodonm9$(EXEEXT): MultiGeneCodonM9.o $(OBJS)
	$(CC) MultiGeneCodonM9.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm8$(EXEEXT): $(PROGSDIR)/multigenecodonm8$(EXEEXT)
$(PROGSDIR)/multigenecodonm8$(EXEEXT): MultiGeneCodonM8.o $(OBJS)
	$(CC) MultiGeneCodonM8.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm2a$(EXEEXT): $(PROGSDIR)/codonm2a$(EXEEXT)
$(PROGSDIR)/codonm2a$(EXEEXT): CodonM2a.o $(OBJS)
	$(CC) CodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcodonm2a$(EXEEXT): $(PROGSDIR)/readcodonm2a$(EXEEXT)
$(PROGSDIR)/readcodonm2a$(EXEEXT): ReadCodonM2a.o $(OBJS)
	$(CC) ReadCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simucodonm2a$(EXEEXT): $(PROGSDIR)/simucodonm2a$(EXEEXT)
$(PROGSDIR)/simucodonm2a$(EXEEXT): SimuCodonM2a.o $(OBJS)
	$(CC) SimuCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm2a$(EXEEXT): $(PROGSDIR)/multigenecodonm2a$(EXEEXT)
$(PROGSDIR)/multigenecodonm2a$(EXEEXT): MultiGeneCodonM2a.o $(OBJS)
	$(CC) MultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

fastreadmultigenecodonm2a$(EXEEXT): $(PROGSDIR)/fastreadmultigenecodonm2a$(EXEEXT)
$(PROGSDIR)/fastreadmultigenecodonm2a$(EXEEXT): FastReadMultiGeneCodonM2a.o $(OBJS)
	$(CC) FastReadMultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenecodonm2a$(EXEEXT): $(PROGSDIR)/readmultigenecodonm2a$(EXEEXT)
$(PROGSDIR)/readmultigenecodonm2a$(EXEEXT): ReadMultiGeneCodonM2a.o $(OBJS)
	$(CC) ReadMultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselm2$(EXEEXT): $(PROGSDIR)/aamutselm2$(EXEEXT)
$(PROGSDIR)/aamutselm2$(EXEEXT): AAMutSelM2.o $(OBJS)
	$(CC) AAMutSelM2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselsparsem9$(EXEEXT): $(PROGSDIR)/aamutselsparsem9$(EXEEXT)
$(PROGSDIR)/aamutselsparsem9$(EXEEXT): AAMutSelSparseM9.o $(OBJS)
	$(CC) AAMutSelSparseM9.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselsparse$(EXEEXT): $(PROGSDIR)/aamutselsparse$(EXEEXT)
$(PROGSDIR)/aamutselsparse$(EXEEXT): AAMutSelSparseOmega.o $(OBJS)
	$(CC) AAMutSelSparseOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselm2$(EXEEXT): $(PROGSDIR)/multigeneaamutselm2$(EXEEXT)
$(PROGSDIR)/multigeneaamutselm2$(EXEEXT): MultiGeneAAMutSelM2.o $(OBJS)
	$(CC) MultiGeneAAMutSelM2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselsparsem9$(EXEEXT): $(PROGSDIR)/multigeneaamutselsparsem9$(EXEEXT)
$(PROGSDIR)/multigeneaamutselsparsem9$(EXEEXT): MultiGeneAAMutSelSparseM9.o $(OBJS)
	$(CC) MultiGeneAAMutSelSparseM9.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselsparse$(EXEEXT): $(PROGSDIR)/multigeneaamutselsparse$(EXEEXT)
$(PROGSDIR)/multigeneaamutselsparse$(EXEEXT): MultiGeneAAMutSelSparseOmega.o $(OBJS)
	$(CC) MultiGeneAAMutSelSparseOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselddp$(EXEEXT): $(PROGSDIR)/aamutselddp$(EXEEXT)
$(PROGSDIR)/aamutselddp$(EXEEXT): AAMutSelDSBDPOmega.o $(OBJS)
	$(CC) AAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutseldp$(EXEEXT): $(PROGSDIR)/aamutseldp$(EXEEXT)
$(PROGSDIR)/aamutseldp$(EXEEXT): AAMutSelSBDPOmega.o $(OBJS)
	$(CC) AAMutSelSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselgdp$(EXEEXT): $(PROGSDIR)/aamutselgdp$(EXEEXT)
$(PROGSDIR)/aamutselgdp$(EXEEXT): AAMutSelGSBDPOmega.o $(OBJS)
	$(CC) AAMutSelGSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

discselac$(EXEEXT): $(PROGSDIR)/discselac$(EXEEXT)
$(PROGSDIR)/discselac$(EXEEXT): DiscSelAC.o $(OBJS)
	$(CC) DiscSelAC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

selempfreq$(EXEEXT): $(PROGSDIR)/selempfreq$(EXEEXT)
$(PROGSDIR)/selempfreq$(EXEEXT): SelACEmpiricalFreq.o $(OBJS)
	$(CC) SelACEmpiricalFreq.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

empselac$(EXEEXT): $(PROGSDIR)/empselac$(EXEEXT)
$(PROGSDIR)/empselac$(EXEEXT): EmpiricalSelACOmega.o $(OBJS)
	$(CC) EmpiricalSelACOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

selac$(EXEEXT): $(PROGSDIR)/selac$(EXEEXT)
$(PROGSDIR)/selac$(EXEEXT): SelACOmega.o $(OBJS)
	$(CC) SelACOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

fastselac$(EXEEXT): $(PROGSDIR)/fastselac$(EXEEXT)
$(PROGSDIR)/fastselac$(EXEEXT): FastSelAC.o $(OBJS)
	$(CC) FastSelAC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

fastempselac$(EXEEXT): $(PROGSDIR)/fastempselac$(EXEEXT)
$(PROGSDIR)/fastempselac$(EXEEXT): FastEmpiricalSelAC.o $(OBJS)
	$(CC) FastEmpiricalSelAC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readaamutselsparse$(EXEEXT): $(PROGSDIR)/readaamutselsparse$(EXEEXT)
$(PROGSDIR)/readaamutselsparse$(EXEEXT): ReadAAMutSelSparseOmega.o $(OBJS)
	$(CC) ReadAAMutSelSparseOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readaamutselddp$(EXEEXT): $(PROGSDIR)/readaamutselddp$(EXEEXT)
$(PROGSDIR)/readaamutselddp$(EXEEXT): ReadAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) ReadAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneaamutselddp$(EXEEXT): $(PROGSDIR)/multigeneaamutselddp$(EXEEXT)
$(PROGSDIR)/multigeneaamutselddp$(EXEEXT): MultiGeneAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) MultiGeneAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediscselac$(EXEEXT): $(PROGSDIR)/multigenediscselac$(EXEEXT)
$(PROGSDIR)/multigenediscselac$(EXEEXT): MultiGeneDiscSelAC.o $(OBJS)
	$(CC) MultiGeneDiscSelAC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneselac$(EXEEXT): $(PROGSDIR)/multigeneselac$(EXEEXT)
$(PROGSDIR)/multigeneselac$(EXEEXT): MultiGeneSelACOmega.o $(OBJS)
	$(CC) MultiGeneSelACOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigeneaamutselddp$(EXEEXT): $(PROGSDIR)/readmultigeneaamutselddp$(EXEEXT)
$(PROGSDIR)/readmultigeneaamutselddp$(EXEEXT): ReadMultiGeneAAMutSelDSBDPOmega.o $(OBJS)
	$(CC) ReadMultiGeneAAMutSelDSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffsel$(EXEEXT): $(PROGSDIR)/diffsel$(EXEEXT)
$(PROGSDIR)/diffsel$(EXEEXT): DiffSel.o $(OBJS)
	$(CC) DiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiffsel$(EXEEXT): $(PROGSDIR)/readdiffsel$(EXEEXT)
$(PROGSDIR)/readdiffsel$(EXEEXT): ReadDiffSel.o $(OBJS)
	$(CC) ReadDiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediffsel$(EXEEXT): $(PROGSDIR)/multigenediffsel$(EXEEXT)
$(PROGSDIR)/multigenediffsel$(EXEEXT): MultiGeneDiffSel.o $(OBJS)
	$(CC) MultiGeneDiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffseldsparse$(EXEEXT): $(PROGSDIR)/diffseldsparse$(EXEEXT)
$(PROGSDIR)/diffseldsparse$(EXEEXT): DiffSelDoublySparse.o $(OBJS)
	$(CC) DiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiffseldsparse$(EXEEXT): $(PROGSDIR)/readdiffseldsparse$(EXEEXT)
$(PROGSDIR)/readdiffseldsparse$(EXEEXT): ReadDiffSelDoublySparse.o $(OBJS)
	$(CC) ReadDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenediffseldsparse$(EXEEXT): $(PROGSDIR)/multigenediffseldsparse$(EXEEXT)
$(PROGSDIR)/multigenediffseldsparse$(EXEEXT): MultiGeneDiffSelDoublySparse.o $(OBJS)
	$(CC) MultiGeneDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenediffseldsparse$(EXEEXT): $(PROGSDIR)/readmultigenediffseldsparse$(EXEEXT)
$(PROGSDIR)/readmultigenediffseldsparse$(EXEEXT): ReadMultiGeneDiffSelDoublySparse.o $(OBJS)
	$(CC) ReadMultiGeneDiffSelDoublySparse.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

siteom$(EXEEXT): $(PROGSDIR)/siteom$(EXEEXT)
$(PROGSDIR)/siteom$(EXEEXT): SiteOmega.o $(OBJS)
	$(CC) SiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenesiteom$(EXEEXT): $(PROGSDIR)/multigenesiteom$(EXEEXT)
$(PROGSDIR)/multigenesiteom$(EXEEXT): MultiGeneSiteOmega.o $(OBJS)
	$(CC) MultiGeneSiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenebranchom$(EXEEXT): $(PROGSDIR)/multigenebranchom$(EXEEXT)
$(PROGSDIR)/multigenebranchom$(EXEEXT): MultiGeneConditionOmega.o $(OBJS)
	$(CC) MultiGeneConditionOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenebranchom$(EXEEXT): $(PROGSDIR)/readmultigenebranchom$(EXEEXT)
$(PROGSDIR)/readmultigenebranchom$(EXEEXT): ReadMultiGeneConditionOmega.o $(OBJS)
	$(CC) ReadMultiGeneConditionOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenesparsebranchom$(EXEEXT): $(PROGSDIR)/multigenesparsebranchom$(EXEEXT)
$(PROGSDIR)/multigenesparsebranchom$(EXEEXT): MultiGeneSparseConditionOmega.o $(OBJS)
	$(CC) MultiGeneSparseConditionOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmultigenesparsebranchom$(EXEEXT): $(PROGSDIR)/readmultigenesparsebranchom$(EXEEXT)
$(PROGSDIR)/readmultigenesparsebranchom$(EXEEXT): ReadMultiGeneSparseConditionOmega.o $(OBJS)
	$(CC) ReadMultiGeneSparseConditionOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ppredtest$(EXEEXT): $(PROGSDIR)/ppredtest$(EXEEXT)
$(PROGSDIR)/ppredtest$(EXEEXT): PostPredTest.o $(OBJS)
	$(CC) PostPredTest.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

.PHONY: format
format:
	clang-format -i *.hpp *.cpp
