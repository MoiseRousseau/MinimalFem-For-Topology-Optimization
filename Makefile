CC=g++
CFLAGS= -Os -I/usr/include/eigen3

all: main.cpp
	$(CC) $(CFLAGS)   main.cpp -o MinimalFEM

clean:
	rm -rf *.o MinimalFEM
