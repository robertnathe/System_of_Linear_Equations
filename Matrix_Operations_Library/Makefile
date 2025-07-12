CXX = g++
CXXFLAGS = -std=c++11 -Wall

all: main

main: main.o MatrixMarket.o
	$(CXX) $(CXXFLAGS) -o main main.o MatrixMarket.o

main.o: main.cpp Matrix.h Vector.h MatrixMarket.h
	$(CXX) $(CXXFLAGS) -c main.cpp

MatrixMarket.o: MatrixMarket.cpp MatrixMarket.h Matrix.h Vector.h
	$(CXX) $(CXXFLAGS) -c MatrixMarket.cpp

clean:
	rm -f *.o main
