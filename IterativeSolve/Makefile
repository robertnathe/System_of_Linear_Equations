IterationSolve: IterationSolve.o
	g++ IterationSolve.o -o IterationSolve
    
IterationSolve.o: IterationSolve.cpp
	g++ -std=c++20 -ggdb -fopenmp -pedantic -c -I -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion IterationSolve.cpp
    
clean:
	rm *.o
