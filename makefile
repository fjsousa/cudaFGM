# The following rules work for UnixWare 2.0.
CC = gcc
NVCC = nvcc

ifdef DEBUG
FLAGS = -g
CUDAFLAGS = -g 
endif

ifdef PROFILING
FLAGS = -g
CUDAFLAGS = -g 
P_FLAG = -pg
endif

ifdef OPTIMIZED
FLAGS = -O3
CUDAFLAGS = -O3
endif

LIBS = -lm
CUDALIBS = -lcudart -lm -lcublas 

%.o : %.c
	$(CC) $(FLAGS) -c $< $(P_FLAG)

%.o : %.cu
	$(NVCC) $(CUDAFLAGS)  -c $< $(P_FLAG)  


OBJ =	fireLib_float.o\
			main.o\
			fireCudaKernel.o


cudaFGM: $(OBJ)
	$(NVCC) -o  $@ $(OBJ) $(CUDALIBS) $(P_FLAG)  

createMaps: createMaps.o 
	$(CC) -o $@ createMaps.o $(LIBS) $(P_FLAG)  

createTestMaps: createTestMaps.o 
	$(CC) -o $@ createTestMaps.o $(LIBS) $(P_FLAG)  

clean:
	@rm -f *.o *.linkinfo *.sw*

cleanDebug:
	@rm -f *.cubin* *.ptx *.gpu* *.cpp* *cudafe* *hash*

FGM_distclean:
	rm cudaFGM  

# End of makefile
