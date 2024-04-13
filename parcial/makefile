# Definición de variables para simplificar cambios futuros
CXX=g++
CXXFLAGS=-O3
parcial_SOURCE=parcial.cpp
parcial_EXECUTABLE=parcial.x

parcialforza_SOURCE=parcialforza.cpp
parcialforza_EXECUTABLE=parcialforza.x


IMAGE_SOURCE=parcialtxt.cpp
IMAGE_EXECUTABLE=parcialtxt.x

.PHONY: all parcial texto clean

# Objetivo por defecto que construye ambos, parcial e 
all: parcial 

# Objetivo para construir y ejecutar la parte de parcial
parcial: $(parcial_SOURCE)
	$(CXX) $(CXXFLAGS) $(parcial_SOURCE) -o $(parcial_EXECUTABLE)
	./$(parcial_EXECUTABLE) | gnuplot

texto: $(IMAGE_SOURCE)
	$(CXX) $(CXXFLAGS) $(IMAGE_SOURCE) -o $(IMAGE_EXECUTABLE)
	./$(IMAGE_EXECUTABLE) > posicion.dat
	gnuplot -p image.gp

parcialforza: $(parcialforza_SOURCE)
	$(CXX) $(CXXFLAGS) $(parcialforza_SOURCE) -o $(parcialforza_EXECUTABLE)
	./$(parcialforza_EXECUTABLE) | gnuplot
# Objetivo para limpiar los archivos generados
clean:
	rm -f $(parcial_EXECUTABLE) $(IMAGE_EXECUTABLE)