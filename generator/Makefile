CFLAGS=-Wall -pedantic -Wno-long-long
CXX=g++
PREP=-c -o
REMOVE=rm -rf
BIN=./generator
DOXYG=Doxyfile
OBJ=./SimConfig.o ./ioproc.o ./generator.o

compile: $(BIN)

clean:
	$(REMOVE) ./*.o $(BIN)

# g++ ./SimConfig.o ./ioproc.o ./generator.o -o ./generator
$(BIN): $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) -o $(BIN)

# g++ SimConfig.cpp -c -o ./SimConfig.o
./SimConfig.o: ./SimConfig.cpp ./SimConfig.h
	$(CXX) $(CFLAGS) ./SimConfig.cpp $(PREP) ./SimConfig.o

# g++ ioproc.cpp -c -o ./ioproc.o
./ioproc.o: ./ioproc.cpp ./ioproc.h ./SimConfig.h
	$(CXX) $(CFLAGS) ./ioproc.cpp $(PREP) ./ioproc.o

# g++ generator.cpp -c -o ./generator.o
./generator.o: ./generator.cpp ./ioproc.h ./SimConfig.h
	$(CXX) $(CFLAGS) ./generator.cpp $(PREP) ./generator.o