# compile command: g++ -o a.out main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11 -fopenmp
CXX=g++
BIN=./simulator
REMOVE=rm -rf
PREP=-c -o
OBJ=./SimConfig.o ./ioproc.o ./main.o

CIMGFLAGS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11
COMPFLAGS=-fopenmp
CFLAGS=-Wall -pedantic -Wno-long-long

compile: $(BIN)

clean:
	$(REMOVE) ./*.o $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) -o $(BIN) $(CIMGFLAGS) $(COMPFLAGS)

./main.o: ./main.cpp
	$(CXX) $(CFLAGS) ./main.cpp $(PREP) ./main.o $(CIMGFLAGS) $(COMPFLAGS)

./SimConfig.o: ./generator/SimConfig.cpp ./generator/SimConfig.h
	$(CXX) $(CFLAGS) ./generator/SimConfig.cpp $(PREP) ./SimConfig.o

./ioproc.o: ./generator/ioproc.cpp ./generator/ioproc.h ./generator/SimConfig.h
	$(CXX) $(CFLAGS) ./generator/ioproc.cpp $(PREP) ./ioproc.o
