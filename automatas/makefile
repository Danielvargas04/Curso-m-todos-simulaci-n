all: karman.gif

forzamiento.x: forzamiento.cpp
	g++ forzamiento.cpp -o forzamiento.x

karman.gif: forzamiento.x
	./forzamiento.x | gnuplot

.PHONY: clean
clean:
	rm -f *.o *~ *.x