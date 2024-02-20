import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

omega = 2.0
dt=0.01

def f1(t, x1, x2):
    return -omega**2*x2

def f2(t, x1, x2):
    return x1

def unpasoRK(t, x1, x2):
    'Realiza un paso de Rungekuta'
    dx11 = dt*f1(t, x1, x2)
    dx12 = dt*f2(t, x1, x2)

    dx21 = dt*f1(t + dt/2, x1 + dx11/2, x2 + dx12/2)
    dx22 = dt*f2(t + dt/2, x1 + dx11/2, x2 + dx12/2)

    dx31 = dt*f1(t + dt/2, x1 + dx21/2, x2 + dx22/2)
    dx32 = dt*f2(t + dt/2, x1 + dx21/2, x2 + dx22/2)

    dx41 = dt*f1(t + dt/2, x1 + dx31, x2 + dx32)
    dx42 = dt*f2(t + dt/2, x1 + dx31, x2 + dx32)
    
    dx1 = (dx11+2*(dx21+dx31)+dx41)/6
    dx2 = (dx12+2*(dx22+dx32)+dx42)/6
        
    return x1 + dx1, x2 + dx2

#Parametros
time = np.arange(0, 7, dt)
x1 = np.zeros(len(time))
x2 = np.zeros(len(time))

#condiciones iniciales
x1[0] = 1
x2[0] = 0

for i in range(len(time)-1):
    x1[i+1], x2[i+1] = unpasoRK(time[i], x1[i], x2[i])

fig, ax = plt.subplots()
ax.plot(time, x2, label = 'Posicion')
ax.plot(time, x1, label = 'Velocidad')
ax.legend()
ax.grid(True)
ax.set_title('Oscilador armonico resuelto con Rungekuta')
plt.show()

