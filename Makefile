CXX = g++
CXXFLAGS = -std=c++17 -g -O0 -Wall -Wextra -pedantic
.PHONY: environment gdb profile clean
environment:
	conda env create -f environment.yml
gdb: gdb_debugging/closest_pair
gdb_debugging/closest_pair: gdb_debugging/closest_pair.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
profile: profiling/distances_before
profiling/distances_before: profiling/profile_distances.cpp
	$(CXX) -std=c++17 -g -O0 -pg -Wall -Wextra $< -o $@
clean:
	rm -f gdb_debugging/closest_pair profiling/distances_before profiling/distances_after profiling/before_release profiling/after_release gmon.out
