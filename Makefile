# compile command: g++ -o a.out main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11 -fopenmp

CC=g++

CFLAGS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11 -fopenmp

build: main

main: main.cpp
	$(CC) -o a.out main.cpp $(CFLAGS)
