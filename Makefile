CC = g++
CXXFLAGS = -O3 -fPIC -fwrapv -fno-strict-aliasing \
	 -Wall -Werror -pedantic -std=c++11 -g -pg \
	-I$(SAGE_LOCAL)/include \
	-I$(SAGE_LOCAL)/include/python2.7
LDFLAGS = -lgsl -lgslcblas -lm -lstdc++

OUTPUT_OPTION = -MMD -MP -o $@
-include $(DEP)

EXEC      = Simulation

SRC = lattice.cpp simulation.cpp explicit_instantiation.cpp
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)

all: phi4.so
	python ./test.py

phi4.cpp: phi4.pyx phi4.pxi
	cython --cplus $<

phi4.o: phi4.cpp

phi4.so: phi4.o $(OBJ)
	$(CC) -shared -pthread -o phi4.so $^ $(LDFLAGS)


test: $(EXEC)
	./Simulation -12600 100 32 32 1000 1000
	./Simulation -12600 100 128 128 1000 1000

$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

clean:
	rm -f phi4.cpp phi4.o phi4.so $(OBJ) $(DEP) $(EXEC)

.PHONY: clean

-include $(DEP)


