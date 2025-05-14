#Build UNS2STR tool, by Mauro Minervino

EXE = UNS2STR
SRC = UNS2STR.cpp
CC = icc
CFLAGS = -std=c++11
TECIOLIBS = /usr/local/apps/tecplot2022r1/360ex_2022r1/bin/libtecio.so

default_target: $(EXE)

.PHONY : clean

clean:
	-rm $(EXE)

$(EXE): UNS2STR.cpp
	$(CC) $(CFLAGS) -o $(EXE) $(SRC) $(TECIOLIBS)
