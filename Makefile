# compile command: g++ -o a.out main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11 -fopenmp
CXX=g++
BIN=./simulator
BIN_SSE=./simulator_sse
BIN_SSE_F=./simulator_sse_fast
REMOVE=rm -rf
PREP=-c -o
OBJ=./SimConfig.o ./ioproc.o ./main.o
OBJ_SSE=./SimConfig.o ./ioproc.o ./main_sse.o
OBJ_SSE_F=./SimConfig.o ./ioproc.o ./main_sse_fast.o

CIMGFLAGS=-O3 -L/usr/X11R6/lib -lm -lpthread -lX11
COMPFLAGS=-fopenmp -fopt-info-vec
CFLAGS=-Wall -pedantic -Wno-long-long

compile: $(BIN) $(BIN_SSE) $(BIN_SSE_F) 

clean:
	$(REMOVE) ./*.o $(BIN) $(BIN_SSE) $(BIN_SSE_F)

$(BIN): $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) -o $(BIN) $(CIMGFLAGS) $(COMPFLAGS)

$(BIN_SSE): $(OBJ_SSE)
	$(CXX) $(CFLAGS) $(OBJ_SSE) -o $(BIN_SSE) $(CIMGFLAGS) $(COMPFLAGS)

$(BIN_SSE_F): $(OBJ_SSE_F)
	$(CXX) $(CFLAGS) $(OBJ_SSE_F) -o $(BIN_SSE_F) $(CIMGFLAGS) $(COMPFLAGS)

./main.o: ./main.cpp
	$(CXX) $(CFLAGS) ./main.cpp $(PREP) ./main.o $(CIMGFLAGS) $(COMPFLAGS)

./main_sse.o: ./main_sse.cpp
	$(CXX) $(CFLAGS) ./main_sse.cpp $(PREP) ./main_sse.o $(CIMGFLAGS) $(COMPFLAGS)

./main_sse_fast.o: ./main_sse_fast.cpp
	$(CXX) $(CFLAGS) ./main_sse_fast.cpp $(PREP) ./main_sse_fast.o $(CIMGFLAGS) $(COMPFLAGS)

./SimConfig.o: ./generator/SimConfig.cpp ./generator/SimConfig.h
	$(CXX) $(CFLAGS) ./generator/SimConfig.cpp $(PREP) ./SimConfig.o

./ioproc.o: ./generator/ioproc.cpp ./generator/ioproc.h ./generator/SimConfig.h
	$(CXX) $(CFLAGS) ./generator/ioproc.cpp $(PREP) ./ioproc.o
