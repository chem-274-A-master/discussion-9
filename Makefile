CXX = g++
CXXFLAGS = -std=c++17 -g -O0 -Wall -Wextra -pedantic
.PHONY: environment gdb clean
environment:
	conda env create -f environment.yml
gdb: gdb_debugging/closest_pair
gdb_debugging/closest_pair: gdb_debugging/closest_pair.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
clean:
	rm -f gdb_debugging/closest_pair
