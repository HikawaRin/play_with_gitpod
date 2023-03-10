CXX = g++
CXXFLAGS = -g -Wall

nsgaii_files = nsgaii.cpp
nsgaii.o:
	-@mkdir build
	cd build/ && $(CXX) $(CXXFLAGS) -c ../$(nsgaii_files) -o nsgaii.o

test_files = test.cpp
exe: nsgaii.o
	-@mkdir build
	cd build/ && $(CXX) $(CXXFLAGS) ../$(test_files) nsgaii.o -o test 

input = in
.PHONY test: exe build
	-@mkdir build
	cd build/ && cp ../$(input) . && ./test <$(input) >out

.PHONY debug: exe build
	cd build/ && cp ../$(input) . && gdb test

.PHONY clean:
	rm -rf build