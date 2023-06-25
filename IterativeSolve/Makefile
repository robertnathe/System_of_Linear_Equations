IterativeSolve: IterativeSolve.o
	g++ IterativeSolve.o -o IterativeSolve
    
IterativeSolve.o: IterativeSolve.cpp
	g++ -std=c++20 -ggdb -fopenmp -pedantic -c -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion IterativeSolve.cpp
    
clean:
	rm *.o
