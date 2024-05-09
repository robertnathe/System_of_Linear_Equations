CC = mpic++
CFLAGS = -Wall -I /usr/lib/x86_64-linux-gnu/
 
# ****************************************************ru
# Targets needed to bring the executable up to date
 
output: main.o
	$(CC) $(CFLAGS) -o output main.o
 
# The main.o target can be written more simply
 
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm *.o

run:
	./output
