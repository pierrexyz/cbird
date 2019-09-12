CXX = icc
CXXFlags = -Wall -Wextra -pedantic -g -O2 -DNDEBUG

EXEC =	cbird
SRC = $(wildcard cnest/*.cpp)
OBJ = $(SRC: .cpp=.o) 
STD = -std=c++0x
LIBPATH = -I/exports/pierre/libraries/include -L/exports/pierre/libraries/lib
GSL = -lm -lgsl -lgslcblas
CUBA = -lcuba -lm
FFTW = -lfftw3


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $^ -o $@ $(STD) $(GSL) $(FFTW) $(LIBPATH) $(CUBA)

%.o: %.cpp
	$(CXX) $(CXXFlags) -o $@ -c $^ 

.PHONY: clean

clean:
	rm $(EXEC) *~

