# trilinos requires boost_system and boost_filesystem
# gnuplotio requires boost_system, boost_filesystem, and boost_iostreams

solver=trilinos
plotter=gnuplot

# boost (note: gnuplot_io depends on boost)
boost_incl = 
boost_lib = -lboost_filesystem -lboost_system -lboost_iostreams -DNDEBUG -DBOOST_UBLAS_NDEBUG

# trilinos
trilinos_incl = -I/usr/include/trilinos
trilinos_lib = -ltrilinos_epetra -ltrilinos_amesos -ltrilinos_teuchoscore

home = $(shell echo ~)

CXX   = g++
FLAGS = -std=c++11 -Wall -O4
LIB   = $(boost_lib)
INC   = -I. $(boost_incl)
DEF   = #-DSOLVER_TYPE=BOOST -DPLOTTER_TYPE=NONE
BIN   = execute

ifeq ($(solver),boost)
	DEF = -DSOLVER_TYPE=BOOST
else ifeq ($(solver),trilinos)
	LIB := $(LIB) $(trilinos_lib)
	INC := $(INC) $(trilinos_incl)
	DEF = -DSOLVER_TYPE=TRILINOS
endif

ifeq ($(plotter),none)
	DEF := $(DEF) -DPLOTTER_TYPE=NONE
else ifeq ($(plotter),gnuplot)
	DEF := $(DEF) -DPLOTTER_TYPE=GNUPLOT
endif

HEADERS = $(shell find . -name "*.h")
SRC = $(shell find . -name "*.cpp")
OBJ = $(patsubst %.cpp, %.o, $(SRC))
DEP = $(patsubst %.cpp, %.d, $(SRC))
TMP = $(shell find . -name "*~")
LOG = $(shell find . -name "*.log")

#--------------------------------------------------------
# Rules for building objects and executable

$(BIN) : $(OBJ) 
	$(CXX) $(FLAGS) -o $(BIN) $(OBJ) $(LIB)

%.o : %.cpp Makefile
	$(CXX) $(FLAGS) $(INC) $(DEF) -c -o $@ $<

#--------------------------------------------------------
# Automatic dependencies: for each %.cpp, make a %.d that
# lists the dependencies of %.cpp, then paste the contents
# of all %.d files into the makefile.

%.d : %.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM -MT $*.o $(FLAGS) $(INC) $< | sed 's|\($*\)\.o[ :]*|\1.o $@ : |g' > $@

-include $(DEP)

#--------------------------------------------------------

trilinos:
	make solver=trilinos

boost:
	make solver=boost

debug:
	make FLAGS='-g -Wall -std=c++0x'

profile:
	make FLAGS='-g -pg -Wall -std=c++0x'

static:
	make FLAGS='-static -Wall -std=c++0x -O4'

cleanall :
	rm -f $(OBJ) $(DEP) $(TMP) $(LOG) $(BIN)

clean :
	rm -f $(OBJ) $(DEP) $(TMP) $(LOG)

cleandep :
	rm -f $(DEP)

cleanobj :
	rm -f $(OBJ)

cleantmp :
	rm -f $(TMP)

new :
	make cleanall
	make


