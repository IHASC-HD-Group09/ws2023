CC = g++

CCFLAGSBASE = -std=c++17 -O1 -funroll-loops -fargument-noalias
CCFLAGS = $(CCFLAGSBASE)

LFLAGS = -lm -lpthread

all: pointer_chasing


pointer_chasing: pointer_chasing.cc Makefile
	$(CC) $(CCFLAGS) -o $@.out $< $(LFLAGS)

clean:
	rm -f pointer_chasing.out