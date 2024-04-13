import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def dispercion(data_path):
    datos=pd.read_csv(data_path, delimiter='\t', names=['t','x'])

    fig, ax = plt.subplots()
    ax.scatter(datos['t'], datos['x'], label = 'Gamma', s=10)
    ax.set_xlabel("Tiempo (t)")
    ax.set_ylabel("Posicion (x)")
    #ax.set_yscale("log")
    ax.legend()
    ax.grid(True)
    ax.set_title('Posicion en funcion del tiempo caso forzado')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python graficador.py <data_path>")
        sys.exit(1)

    data_path = sys.argv[1]
    
    dispercion(data_path)