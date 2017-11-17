import numpy as np
import matplotlib.pyplot as plt
import geral.runge

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


