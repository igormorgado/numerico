import math
import bitstring
import matplotlib as mpl
import matplotlib.pyplot as plt
from erros import eps


def err_PF(n, x):
    """Propaga o erro de ponto flutuante para n iteracoes sobre x"""
    soma = [x]
    prod = [x]
    result = []

    for i in range(1, n):
        soma.append(soma[-1] + x)
        prod.append((i + 1) * x)

    for i in range(n):
        result.append(abs(soma[i] - prod[i]))

    return result


def err_plot(n, x, title=None):
    """Gera o gráfico de erros para uma lista de valores iniciais"""
    if isinstance(x, list) or isinstance(x, tuple):
        result = []
        for i in x:
            lab = str(i)
            if len(lab) > 7:
                lab = lab[:7] + "..."
            result.append(plt.plot(range(n), err_PF(n, i), alpha=0.5, label=lab))
    else:
        lab = str(i)
        if len(lab) > 7:
            lab = lab[:7] + "..."

        result = plt.plot(range(n), err_PF(n, x), alpha=0.5, label=lab)

    if title is not None:
        plt.title(title)

    plt.legend()

    return result


def erro_local(l):
    """Retorna a lista com o gradiente por iteracao"""
    result = []
    for i in range(1, len(l)):
        result.append(abs(l[i] - l[i - 1]))

    return result


def erro_global(l):
    """Retorna a lista com o gradiente por iteracao"""
    result = []
    for i in range(1, len(l)):
        result.append(abs(l[i] - l[0]))

    return result


def err_step(n, x):
    """Gera o gráfico dos erros para cada iteração
    """


# 'double' means IEEE 754 double precision -- c 'double'
epsilon = math.ldexp(1.0, -53)  # smallest double that 0.5+epsilon != 0.5
maxDouble = float(2**1024 - 2**971)  # From the IEEE 754 standard
minDouble = math.ldexp(1.0, -1022)  # min positive normalized double
smallEpsilon = math.ldexp(1.0, -1074)  # smallest increment for doubles < minFloat
infinity = math.ldexp(1.0, 1023) * 2


def nextafter(x, y):
    """returns the next IEEE double after x in the direction of y if possible"""
    if y == x:
        return y  # if x==y, no increment

    # handle NaN
    if x != x or y != y:
        return x + y

    if x >= infinity:
        return infinity

    if x <= -infinity:
        return -infinity

    if -minDouble < x < minDouble:
        if y > x:
            return x + smallEpsilon
        else:
            return x - smallEpsilon

    m, e = math.frexp(x)
    if y > x:
        m += epsilon
    else:
        m -= epsilon

    return math.ldexp(m, e)


def reprbin(x):
    """Retorna a representação strem de bit na forma (sinal, expoente, mantissa)"""
    b = bitstring.pack('>d', x)
    s, w, p = b[:1], b[1:12], b[12:]
    return s, w, p
