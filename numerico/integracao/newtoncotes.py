import numpy as np
from numerico.interpolacao.lagrange import interpola


# Integracao newton-cotes

def retangulos(a, b, x, y, n=2):
    """Integra o intervalo [a,b] por retangulos  aproximando-se nos pontos 
    (x,y) utilizando n subintervalos"""

    t, h = np.linspace(a, b, n, retstep=True)
    f_t =  interpola(t, x, y)

    soma = 0
    for i in range(1,len(f_t)-1):
        soma += f_t[i] * h 
    
    return soma


def trapezio(a, b, x, y, n=2):
    """Integra o intervalo [a,b] por trapezios  aproximando-se nos pontos 
    (x,y) utilizando n subintervalos"""

    t, h = np.linspace(a, b, n, endpoint=True, retstep=True)
    f_t =  interpola(t, x, y)

    soma = 0
    for i in range(len(f_t)-1):
        soma += (f_t[i]+f_t[i+1]) * h / 2
    
    # sum((f_t[1:] + f_t[:-1]) * (h/2))
    return soma

def simpson(a, b, x, y, n=2):
    """Integra o intervalo [a,b] por polinomio de segunda ordem aproximando-se
    nos pontos (x,y) utilizando n subintervalos"""

    if ( n % 2 != 0 ):
        raise ValueError("Numero de pontos deve ser par")

    t, h = np.linspace(a, b, n, endpoint=True, retstep=True)
    f_t =  interpola(t, x, y)

    soma = 0
    for i in range(0,len(f_t)-2,2):
        soma += (f_t[i] + 4*f_t[i+1] + f_t[i+2]) * h / 3
    
    # sum((f_t[1:] + f_t[:-1]) * (h/2))
    return soma

