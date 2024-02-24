'''Programa pensado para graficar los diferentes datos de los programas en c++ 
para no usar gnuplot si no las ventajas de python y matplotlib'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

datos=pd.read_csv('data/datosparabola.dat',sep=' ', names=['x','y'])
#datoskuta=pd.read_csv('data/datoskuta.dat',sep=' ', names=['x','y','z'])

fig, ax = plt.subplots()
ax.plot(datos['x'], datos['y'], label = 'Posicion')
ax.legend()
ax.grid(True)
ax.set_title('Oscilador armonico resuelto con Rungekuta')
plt.show()
