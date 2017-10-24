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
CPPFLAGS= -std=c++0x -g -Wall -O3 $(USE_GSL) $(INCLUDES)
LDFLAGS=$(STATIC)
INSTALL_DIR=
INSTALL_LIB=
SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp AAMutSelOmegaCodonSubMatrix.cpp GTRSubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp linalg.cpp cdf.cpp Chain.cpp MultiGeneChain.cpp MPIBuffer.cpp MultiGeneMPIModule.cpp CodonM2aModel.cpp MultiGeneCodonM2aModel.cpp AAMutSelM2aModel.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= aamutselm2a aamutseldp dpom dgamom gammasiteom diffsel multigenecodonm2a codonm2a aamutsel multigeneglobom globom chi2
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

chi2$(EXEEXT): $(PROGSDIR)/chi2$(EXEEXT)
$(PROGSDIR)/chi2$(EXEEXT): Chi2.o $(OBJS)
	$(CC) Chi2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutseldp$(EXEEXT): $(PROGSDIR)/aamutseldp$(EXEEXT)
$(PROGSDIR)/aamutseldp$(EXEEXT): AAMutSelSBDPOmega.o $(OBJS)
	$(CC) AAMutSelSBDPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutsel$(EXEEXT): $(PROGSDIR)/aamutsel$(EXEEXT)
$(PROGSDIR)/aamutsel$(EXEEXT): AAMutSelOmega.o $(OBJS)
	$(CC) AAMutSelOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aamutselm2a$(EXEEXT): $(PROGSDIR)/aamuselm2a$(EXEEXT)
$(PROGSDIR)/aamutselm2a$(EXEEXT): AAMutSelM2a.o $(OBJS)
	$(CC) AAMutSelM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm2a$(EXEEXT): $(PROGSDIR)/codonm2a$(EXEEXT)
$(PROGSDIR)/codonm2a$(EXEEXT): CodonM2a.o $(OBJS)
	$(CC) CodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonm8$(EXEEXT): $(PROGSDIR)/codonm8$(EXEEXT)
$(PROGSDIR)/codonm8$(EXEEXT): CodonM8.o $(OBJS)
	$(CC) CodonM8.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm2a$(EXEEXT): $(PROGSDIR)/multigenecodonm2a$(EXEEXT)
$(PROGSDIR)/multigenecodonm2a$(EXEEXT): MultiGeneCodonM2a.o $(OBJS)
	$(CC) MultiGeneCodonM2a.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigenecodonm8$(EXEEXT): $(PROGSDIR)/multigenecodonm8$(EXEEXT)
$(PROGSDIR)/multigenecodonm8$(EXEEXT): MultiGeneCodonM8.o $(OBJS)
	$(CC) MultiGeneCodonM8.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

multigeneglobom$(EXEEXT): $(PROGSDIR)/multigeneglobom$(EXEEXT)
$(PROGSDIR)/multigeneglobom$(EXEEXT): MultiGeneSingleOmega.o $(OBJS)
	$(CC) MultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

dgamom$(EXEEXT): $(PROGSDIR)/dgamom$(EXEEXT)
$(PROGSDIR)/dgamom$(EXEEXT): DiscGammaSiteOmega.o $(OBJS)
	$(CC) DiscGammaSiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

dpom$(EXEEXT): $(PROGSDIR)/dpom$(EXEEXT)
$(PROGSDIR)/dpom$(EXEEXT): DPOmega.o $(OBJS)
	$(CC) DPOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

globom$(EXEEXT): $(PROGSDIR)/globom$(EXEEXT)
$(PROGSDIR)/globom$(EXEEXT): SingleOmega.o $(OBJS)
	$(CC) SingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

gammasiteom$(EXEEXT): $(PROGSDIR)/gammasiteom$(EXEEXT)
$(PROGSDIR)/gammasiteom$(EXEEXT): GammaSiteOmega.o $(OBJS)
	$(CC) GammaSiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffsel$(EXEEXT): $(PROGSDIR)/diffsel$(EXEEXT)
$(PROGSDIR)/diffsel$(EXEEXT): DiffSel.o $(OBJS)
	$(CC) DiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

