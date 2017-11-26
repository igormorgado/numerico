import numpy as np
from numerico.diferencas.dif_divididas import dif_div_prog


def coeficiente_interpolacao(x,y):
    """Returns the highest order interpolation coefficients 
    to a set of data points"""

    return dif_div_progr(x,y)[:,0]

def interpola(t, x, y):
    """Returns the points 't' interpolated by 'x,y' points.
    The interpolation order will be based on size of 'x,y'"""

    b = coeficiente_interpolacao(x, y)

    f = np.full_like(t, b[0], dtype='float64')

    for i in range(1,len(b)):
        mul = 1
        for m in range(i):
            mul *= t - x[m]

        f += b[i] * mul

    return f
