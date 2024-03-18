
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /usr/lib/x86_64-linux-gnu/hdf5/serial
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/homebrew/opt/hdf5
endif

CC = clang
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib

OBJ = main.o readpar.o readdata.o makemap.o

default: cart2pol

%.o: %.c defs.h
	$(CC) $(FLAGS) $(INC) -c $<

cart2pol: $(OBJ) defs.h
	$(CC) $(FLAGS) $(LIB) -o cart2pol $(OBJ) -lhdf5 -lm

clean:
	rm -f *.o cart2pol

cleanall:
	rm -f *.o *.out cart2pol
