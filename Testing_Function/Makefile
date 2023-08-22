# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -std=c++20 -ggdb -fopenmp -pedantic -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion
 
# ****************************************************ru
# Targets needed to bring the executable up to date
 
output: main.o GeneralPurpose.o ArithmeticRoutines.o
	$(CC) $(CFLAGS) -o output main.o GeneralPurpose.o ArithmeticRoutines.o
 
# The main.o target can be written more simply
 
main.o: main.cpp GeneralPurpose.h ArithmeticRoutines.h
	$(CC) $(CFLAGS) -c main.cpp
 
GeneralPurpose.o: GeneralPurpose.h
 
ArithmeticRoutines.o: ArithmeticRoutines.h GeneralPurpose.h

clean:
	rm *.o

run: 
	./output
