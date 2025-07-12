CC = g++
CFLAGS = -std=c++20 -ggdb -fopenmp -pedantic -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

output: main.o
	$(CC) $(CFLAGS) -o output main.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm *.o output

run:
	./output
