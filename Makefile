CXX      ?= g++
CXXFLAGS ?= -O2 -fopenmp
SRC      := $(wildcard src/*.cpp)

secretscores: $(SRC) $(wildcard src/*.h)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) -lz

clean:
	rm -f secretscores
