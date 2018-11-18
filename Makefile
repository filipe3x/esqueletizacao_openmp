# source files.
SRC = img.cpp skeletonize.cpp ppm.cpp utils.cpp papi_inst.cpp main.cpp

OBJ = $(SRC:.cpp=.o)

OUT = skeletonize

# include directories
INCLUDES = -I. -I/usr/local/papi/include
 
# C++ compiler flags (-g -O2 -Wall)
#CCFLAGS = -O2 -Wall -static
#CCFLAGS = -O3 -fopenmp -mavx -march=native
#CCFLAGS = -O0
CCFLAGS = -O3 -march=native -fopenmp
#CCFLAGS = -O3

# compiler
CCC = g++-6 
#CCC = /opt/intel/Compiler/11.1/073/bin/ia32/icpc 
#CCC = g++-4.5
# library paths
LIBS = -L/usr/local/papi/lib -lm -lpapi

# compile flags
LDFLAGS = -g

.SUFFIXES: .cpp .c 


default: $(OUT)

.cpp.o:
	$(CCC) $(CCFLAGS) $(INCLUDES)  -c $< -o $@

.c.o:
	$(CCC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

$(OUT): $(OBJ)
	$(CCC) -o $(OUT) $(CCFLAGS) $(OBJ) $(LIBS) 

depend:  dep
#
#dep:
#	makedepend -- $(CFLAGS) -- $(INCLUDES) $(SRC)

clean:
	rm -f *.o .a *~ Makefile.bak 
