all: karman.gif

forzamiento.x: forzamiento.cpp
	g++ forzamiento.cpp -o forzamiento.x

karman.gif: forzamiento.x
	./forzamiento.x | gnuplot

ondas:
	g++ luz_ondas.cpp -o luz_ondas.x
	./luz_ondas.x | gnuplot
.PHONY: clean
clean:
	rm -f *.o *~ *.x