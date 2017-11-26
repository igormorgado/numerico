import numpy as np

def vandermonde(x):
    """Retorna uma matriz de vandermonde dada pelos pontos 'x'"""
    
    n = len(x)
    v = np.tile(x,(n,1))
    e = np.arange(n)
    v = v.T ** e

    return v

