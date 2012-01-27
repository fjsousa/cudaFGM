# The following rules work for UnixWare 2.0.
CC = gcc
NVCC = nvcc

ifdef DEBUG
FLAGS = -g
CUDAFLAGS = -g
endif

ifdef OPTIMIZED
FLAGS = -O3
CUDAFLAGS = -O3
endif

LIBS = -lm
CUDALIBS = -lcudart -lm -lcublas 

%.o : %.c
	$(CC) $(FLAGS) -c $< 

%.o : %.cu
	$(NVCC) $(CUDAFLAGS)  -c $< 


OBJ =	fireLib_float.o\
			main.o\
			FK_SpreadAtNeighbors.o


cudaFGM: $(OBJ)
	$(NVCC) -o  $@ $(OBJ) $(CUDALIBS) 

createMaps: createMaps.o 
	$(CC) -o $@ createMaps.o $(LIBS) 

createTestMaps: createTestMaps.o 
	$(CC) -o $@ createTestMaps.o $(LIBS) 

clean:
	@rm -f *.o *.linkinfo *.sw*

cleanDebug:
	@rm -f *.cubin* *.ptx *.gpu* *.cpp* *cudafe* *hash*

FGM_distclean:
	rm cudaFGM  

# End of makefile
