CXX = g++
CXXFlags = -Wall -Wextra -pedantic -g -O2 -DNDEBUG

EXEC =	cbird
SRC = $(wildcard *.cpp)
OBJ = $(SRC: .cpp=.o) 
STD = -std=c++11
GSL = -lm -lgsl -lgslcblas
FFTW = -lfftw3


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $^ -o $@ $(STD) $(GSL) $(FFTW)
	
%.o: %.cpp
	$(CXX) $(CXXFlags) -o $@ -c $^ 

.PHONY: clean

clean:
	rm $(EXEC) *~

