import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


def spline_cubico_coef(x,y):
    """Retorna uma lista com os coeficientes do spline cubico.
    A entrada são as listas com o conjunto de pontos (x, f(x))
    conhecidos da uma funcao dada."""
    
    N = len(x) - 1
    
    # Calculando as distancias dos pontos
    h = x[1:] - x[:-1]

    ##########################
    # Primeira derivada
    yp = (y[1:] - y[:-1])/h

    ##########################
    # Diferenca das derivadas
    dyp = np.zeros(N+1)
    dyp[1:-1] = (yp[1:] - yp[:-1])

    # Determinando os coeficientes de 'a'
    a = y.copy()
    
    # Matriz de Distancias

    # Diagonal Principal
    dp = np.ones(N+1)
    dp[1:-1] = 2 * (h[1:] + h[:-1])

    # Diagonal superior
    du = np.zeros(N)
    du[1:] = h[1:]

    # Diagonal inferior
    dl = np.zeros(N)
    dl[:-1] = h[:-1]

    H =  np.diag(dp)
    H += np.diag(du,k=1)
    H += np.diag(dl,k=-1)

    # Resolvendo os coeficientes C 
    c = la.solve(H,3*dyp)
    
    # Resolvendo os coeficientes D
    d = (c[1:] - c[:-1])/(3*h)
    
    # Resolvendo os coeficientes B
    b = yp - c[:-1]*h - d*h**2
    
    return a[:-1], b, c[:-1], d, x


def spline_cubico(x, coef):
    """Retorna of valores interpolados do vetor 'x' baseados em coef"""

    a, b, c, d, t  = coef

    pi = np.array([])
    
    # Valores menores
    s = lambda x: a[0] + b[0]*(x-t[0]) + c[0]*(x-t[0])**2 + d[0]*(x-t[0])**3
    x_ = x[x < t[0]]
    y_ = s(x_)
    pi = np.concatenate((pi,y_))

    for i in range(len(t)-1):
        # Define o polinomio baseado nos coeficientes 
        s = lambda x: a[i] + b[i]*(x-t[i]) + c[i]*(x-t[i])**2 + d[i]*(x-t[i])**3

        # Calcula os valores baseado no intervalo
        x_ = x[(t[i] <= x) & ( x < t[i+1])]
        y_ = s(x_)
        pi = np.concatenate((pi,y_))

    # Valores maiores
    x_ = x[x >= t[-1]]
    y_ = s(x_)
    pi = np.concatenate((pi,y_))
        
    return pi


if __name__ == "__main__":
    x_20 = np.linspace(-5, 5, 20)
    x = x_20
    y = runge(x)

    # Calcula os coeficientes do Spline Cubico
    coef = spline_cubico_coef(x, y)


    a, b, c, d, t  = coef

    for i in range(len(coef[1])):
        polystr = "S[{}] = {:7.4f} + {:7.4f}x + {:7.4f}x^2 + {:7.4f}x^3"
        print("Intervalo {}: [ {}, {} ]".format(i, t[i], t[i+1]))
        print(polystr.format(i, a[i],b[i],c[i],d[i]))
        print()

    X = np.linspace(-5,5, 31)
    Y = spline_cubico(X, coef[0])


    rx = np.linspace(-5,5, 101)

    # Plota a imagem
    plt.plot(x, y, marker='o', linestyle='', color='r', label='Pontos')
    plt.plot(rx, runge(rx), color='k', label='Funcao Real')
    plt.plot(X, Y, color='b', label='Spline')
    plt.ylim(0,1)
    plt.legend()
    plt.show()


