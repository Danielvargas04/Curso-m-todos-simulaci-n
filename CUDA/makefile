#variables
CXX= g++
CXXFLAGS= -O3

CUDA=nvcc
CUDAFLAGS= -arch=sm_75 -gencode=arch=compute_75,code=sm_75

#Objetivo principal
all: main

#Objetivos secundarios
cuda_target=ondas.xu

# Compila los archivos fuente en objetos
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@


%.xu: %.cu 
	$(CUDA) $(CUDAFLAGS) $< -o $(cuda_target)

run-ondas: ondas.xu
	./$(cuda_target)
clean:
	rm -f *.o *.xu main $(cuda_target)

# Pseudotarget para evitar conflictos con archivos del mismo nombre
.PHONY: all clean