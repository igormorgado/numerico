import numpy as np
from numerico.matriz.vandermonde import vandermonde


def coeficientes(x, y):
    """Resolve o sistema linear para encontrar os coeficientes
    interpoladores"""

    return np.linalg.solve(vandermonde(x), y)

def interpola(x, A):
    """Interpola os pontos 'x' com os coeficientes 'A' em um polinomio
    crescente em grau
    
    p_n(x) = a_0 x^0 + a_1 x^1 + a_2 x^2 + ... + a_n x^n
    """

    p = np.zeros_like(x)

    for n,a in enumerate(A):
        p += a * x**n

    return p

