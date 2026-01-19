#Build UNS2STR tool, by Mauro Minervino

EXE = UNS2STR
SRC = UNS2STR.cpp
CC = icpx
CFLAGS = -std=c++11
TECIOLIBS = /opt/tecplot/360ex_2025r2/bin/libtecio.so

default_target: $(EXE)

.PHONY : clean

clean:
	-rm $(EXE)

$(EXE): UNS2STR.cpp
	$(CC) $(CFLAGS) -o $(EXE) $(SRC) $(TECIOLIBS)
