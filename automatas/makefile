all: karman.gif

fluidos.x: fluidos.cpp
	g++ fluidos.cpp -o fluidos.x

karman.gif: fluidos.x
	./fluidos.x | gnuplot

.PHONY: clean
clean:
	rm -f *.o *~ *.x