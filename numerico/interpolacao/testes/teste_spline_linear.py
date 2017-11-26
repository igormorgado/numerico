import numpy as np
from interpolacao.spline_linear import spline_linear
from geral.runge import runge, plota_runge

# Espaco de testes
x_05 = np.linspace(-5,5,5)
x_10 = np.linspace(-5,5,10)
x_20 = np.linspace(-5,5,20)

plota_runge(x_05, legenda='5 pontos')
plota_runge(x_10, legenda='10 pontos')
plota_runge(x_20, legenda='20 pontos')
plota_spline_linear(x_05, runge, 100, legenda='Interpolacao Linear')
plt.legend()
plt.show()


