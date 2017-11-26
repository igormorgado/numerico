import numpy as np

def peso(t, x, j):
    """Retorna o j-esimo peso de Lagrange para os 't' pontos proximos ao
    conjunto x"""

    l = 1
    k = len(x)
    for i in range(k):
        if i != j: 
            l *= (t-x[i])/(x[j] - x[i])

    return l

def interpola(t, x, y):
    """Retorna os pontos 't' interpolados pelo Polinomio de lagrange"""

    L = np.zeros_like(t)
    k = len(x)
    l = peso

    for j in range(k):
        L += y[j] * l(t, x, j)
    
    return L



