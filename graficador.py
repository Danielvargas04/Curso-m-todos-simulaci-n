'''Programa pensado para graficar los diferentes datos de los programas en c++ 
para no usar gnuplot si no las ventajas de python y matplotlib'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

datos=pd.read_csv('data/orbitas.dat',sep=' ', names=['x1','y1', 'x2','y2'])
#datoskuta=pd.read_csv('data/datoskuta.dat',sep=' ', names=['x','y','z'])

fig, ax = plt.subplots()
ax.plot(datos['x1'], datos['y1'], label = 'Planeta 1')
ax.plot(datos['x2'], datos['y2'], label = 'Planeta 2')
ax.legend()
ax.grid(True)
ax.set_title(' ')
plt.show()
