# source files.
SRC = img.cpp skeletonize.cpp ppm.cpp utils.cpp papi_inst.cpp mpi_inst.cpp main.cpp

OBJ = $(SRC:.cpp=.o)

OUT = skeletonize

PAPI = 5.4.1

#DEFINEMACRO = -D PAPI -D TESTING
#DEFINEMACRO = -D PAPI -D PRODUCTION
DEFINEMACRO = -D TESTING -D STATIC

# include directories
#INCLUDES = -I. -I/share/apps/papi/$(PAPI)/include
#INCLUDES = -I. -I/usr/local/papi/include
INCLUDES = -I. -I/usr/local/include/
#INCLUDES = -I. -I/usr/include
 
# C++ compiler flags (-g -O2 -Wall)
#CCFLAGS = -O2 -Wall -static
#CCFLAGS = -O3 -fopenmp -mavx -march=native
#CCFLAGS = -O0
#CCFLAGS = -g -pthread -Wall -march=native -fopenmp
CCFLAGS = -O3 -Wall -march=native --std=c++11 -fopenmp
#CCFLAGS = -Wall -march=native --std=c++11 -ftree-vectorizer-verbose=3 -ftree-vectorize -ftree-slp-vectorize -fopenmp -fopenmp-simd
#CCFLAGS = -O3
#CCFLAGS = -O3 -S -Wall -march=native -fopenmp -fopenmp-simd

# compiler
#CCC = g++-7
CCC = mpic++
#CCC = /opt/intel/Compiler/11.1/073/bin/ia32/icpc 
#CCC = g++-4.5

# library paths
#LIBS = -L/share/apps/papi/$(PAPI)/lib -lm -lpapi
LIBS = -L/usr/local/lib/ -lm -lpapi
#LIBS = -L/usr/local/papi/lib -lm -lpapi
#LIBS = -L/usr/lib/x86_64-linux-gnu -lm -lpapi

# compile flags
LDFLAGS = -g

.SUFFIXES: .cpp .c 


default: $(OUT)

.cpp.o:
	$(CCC) $(DEFINEMACRO) $(CCFLAGS) $(INCLUDES)  -c $< -o $@

.c.o:
	$(CCC) $(DEFINEMACRO) $(CCFLAGS) $(INCLUDES) -c $< -o $@

$(OUT): $(OBJ)
	$(CCC) -o $(OUT) $(CCFLAGS) $(OBJ) $(LIBS) 

depend:  dep
#
#dep:
#	makedepend -- $(CFLAGS) -- $(INCLUDES) $(SRC)

dynamic: DEFINEMACRO = -D TESTING -D DYNAMIC
dynamic: $(OUT)

papi: DEFINEMACRO += -D PAPI
papi: $(OUT)

prod: DEFINEMACRO = -D PRODUCTION -D $(schedule)
prod: $(OUT)

clean:
	rm -f *.o .a *~ Makefile.bak 
