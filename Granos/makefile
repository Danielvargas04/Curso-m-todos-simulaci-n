# Definición de variables para simplificar cambios futuros
CXX=g++
CXXFLAGS=-O3
granos_SOURCE=parcial.cpp
granos_EXECUTABLE=parcial.x
IMAGE_SOURCE=granos_image.cpp
IMAGE_EXECUTABLE=granos_image.x

.PHONY: all granos  clean

# Objetivo por defecto que construye ambos, granos e 
all: granos 

# Objetivo para construir y ejecutar la parte de granos
granos: $(granos_SOURCE)
	$(CXX) $(CXXFLAGS) $(granos_SOURCE) -o $(granos_EXECUTABLE)
	./$(granos_EXECUTABLE) | gnuplot
	./$(granos_EXECUTABLE) > script.gp

# Objetivo para limpiar los archivos generados
clean:
	rm -f $(granos_EXECUTABLE) $(IMAGE_EXECUTABLE)