
RESTART   = h5in
BOUNDARY  = dirichlet#cornermirror
INITIAL   = sedov2d
HYDRO     = euler
RIEMANN   = hllc_lm
GRAVITY   = slab
NOZZLE    = wind
SPACE     = plm
TIMESTEP  = rk2
MESHSPEED = broken_pl#homologous
OUTPUT    = h5out#_withghost

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /usr/lib/x86_64-linux-gnu/hdf5/serial
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/homebrew/opt/hdf5
endif

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib

OBJ = main.o domain.o gridsetup.o mpisetup.o readpar.o exchange.o onestep.o report.o profiler.o $(RESTART).o $(BOUNDARY).o $(INITIAL).o $(HYDRO).o $(RIEMANN).o $(GRAVITY).o $(NOZZLE).o $(SPACE).o $(TIMESTEP).o $(OUTPUT).o $(MESHSPEED).o

default: cube

%.o: %.c defs.h
	$(CC) $(FLAGS) $(INC) -c $<

$(RESTART).o : Restart/$(RESTART).c defs.h
	$(CC) $(FLAGS) $(INC) -c Restart/$(RESTART).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c defs.h
	$(CC) $(FLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(INITIAL).o : Initial/$(INITIAL).c defs.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c defs.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(RIEMANN).o : Riemann/$(RIEMANN).c defs.h
	$(CC) $(FLAGS) $(INC) -c Riemann/$(RIEMANN).c

$(GRAVITY).o : Gravity/$(GRAVITY).c defs.h
	$(CC) $(FLAGS) $(INC) -c Gravity/$(GRAVITY).c

$(NOZZLE).o : Nozzle/$(NOZZLE).c defs.h
	$(CC) $(FLAGS) $(INC) -c Nozzle/$(NOZZLE).c

$(SPACE).o : Space/$(SPACE).c defs.h
	$(CC) $(FLAGS) $(INC) -c Space/$(SPACE).c

$(TIMESTEP).o : Timestep/$(TIMESTEP).c defs.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(MESHSPEED).o : Meshspeed/$(MESHSPEED).c defs.h
	$(CC) $(FLAGS) $(INC) -c Meshspeed/$(MESHSPEED).c

$(OUTPUT).o : Output/$(OUTPUT).c defs.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

cube: $(OBJ) defs.h
	$(CC) $(FLAGS) $(LIB) -o cube $(OBJ) -lhdf5 -lm

clean:
	rm -f *.o cube

rmdata:
	rm -f *.h5

cleanall:
	rm -f *.o *.h5 *.dat cube
