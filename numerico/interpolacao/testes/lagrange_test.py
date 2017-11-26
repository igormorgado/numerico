from numerico.interpolacao.lagrange import *
import matplotlib.pyplot as plt
import numpy as np


# Parabola
#x = np.array([0,10,15,22.5,30])
#y = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.36])

x = np.array([0,10,15,22.5,30])
y = np.array([1.1,3.8,2.1,7.4, -3.2])

t = np.linspace(x[0], x[-1], 101) 

p = interpola(t, x, y)

plt.plot(x, y, 'ro', label="orig")
plt.plot(t, p, 'b-', label="interp")
plt.legend()
plt.show()

