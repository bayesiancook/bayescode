CC=$(CXX)
CPPFLAGS= -std=c++11 -Wall -Wextra -Wno-unused-parameter -O3
# CPPFLAGS= -std=c++11 -g -Wall -Wextra -Wno-unused-parameter -fprofile-arcs -ftest-coverage
LDFLAGS =
# LDFLAGS = -fprofile-arcs -ftest-coverage
SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp GTRSubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp Chain.cpp AAMutSelOmegaCodonSubMatrix.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= globom
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean format
all: $(PROGS)

# Rules to generate the dependencies automatically
%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$

# Rules to create .o files from .cpp files
%.o: %.cpp %.d tinycompo.hpp
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets
$(PROGSDIR)/globom: SingleOmega.o $(OBJS)

$(PROGS):
	$(CC) $^ $(LDFLAGS) -o $@

clean:
	-rm -f *.o *.d *.d.* *.gcno codonm8/*.o codonm8/*.d *.gcda *.gcov
	-rm -f $(PROGS)
	-rm -rf html/ latex/
	-rm -f tinycompo.hpp

format:
	clang-format -i $(ALL_SRCS) *.hpp
