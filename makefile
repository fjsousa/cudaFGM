# makefile for fireSim, a sample fire behavior simulator using fireLib
# Collin D. Bevins, October 1996


#ifdef $(DEVEMULATION)
# DEV_EMULATION=-deveemu
#endif

# The following rules work for UnixWare 2.0.
#CC = gcc
NVCC = nvcc
CUDAFLAGS = -arch=sm_13 -O3 

LIBS = -lm
CUDALIBS = -lcudart -lm#-lcublas -lcudart 

%.o : %.cu
	$(NVCC) $(CUDAFLAGS)  -c $< 


OBJ =	fireCudaLib.o\
			fireWrapper.o\
			FK_NoWindNoSlope.o\
			FK_SpreadAtNeighbors.o\
			FK_WindAndSlope.o\
			errorStuff.o\
			createMaps.o


fireCuda: $(OBJ)
	$(NVCC) -o  $@ $(OBJ) $(CUDALIBS) 

clean:
	@rm -f *.o *.linkinfo *.sw*

cleanDebug:
	@rm -f *.cubin* *.ptx *.gpu* *.cpp* *cudafe* *hash*

distclean: clean
	rm fireSim

cleandist: distclean


# End of makefile
