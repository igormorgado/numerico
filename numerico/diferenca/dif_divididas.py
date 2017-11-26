import numpy as np

def dif_prog(y):
    """Retorna a diferenca progressiva de um numero de pontos"""

    dd = np.zeros([len(y), len(y)]) 
    dd[0] = y
    dd[1,:-1] = (dd[0,1:] - dd[0,:-1])

    for i in range(2,len(y)):
        dd[i,:-i] = dd[i-1,1:-(i-1)] - dd[i-1,:-i]

    return dd


def dif_prog_div(x, y):
    """Retorna a diferenca progressiva dividida de um numero de pontos"""

    dd = np.zeros([len(y), len(y)]) 
    dd[0] = y
    dd[1,:-1] = (dd[0,1:] - dd[0,:-1]) / (x[1:]-x[:-1])

    for i in range(2,len(y)):
        dd[i,:-i] = (dd[i-1,1:-(i-1)] - dd[i-1,:-i]) / (x[i:]-x[:-i])

    return dd


