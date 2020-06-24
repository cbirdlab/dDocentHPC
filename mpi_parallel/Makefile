SRCS = mpi_parallel.cpp
OBJS = $(SRCS:.cpp=.o)
BIN  = mpi_parallel

# openmpi mpich ...
# CXX =  mpic++
# LD  =  mpic++

# intel mpi
CXX =  mpiicpc
LD  =  mpiicpc

CXXFLAGS = -std=c++11 -Wall
LDFLAGS  = 

all: $(BIN)


$(BIN): $(OBJS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(BIN)

clean:
	rm $(OBJS) $(BIN)

.o:
	$(CXX) $(CXXFLAGS) -c $<

