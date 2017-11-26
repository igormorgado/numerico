#!/usr/bin/env python3

import math
from scipy import optimize
from sympy.abc import x
from sympy import diff
from sympy import lambdify

from bisseccao import bisseccao
from falsaposicao import falsaposicao
from pontofixo import pontofixo1
from pontofixo import pontofixo2
from newtonraphson import newtonraphson
from secante import secante

eps = 1e-15


def f1(x):
    """x * log(x) - 1"""
    return x * math.log(x) - 1


def f2(x):
    """x**3 - 9*x + 3"""
    return x**3 - 9 * x + 3


def phi(x):
    """x**3/9 + 1/3"""
    return x**3 / 9 + 1 / 3


def f3(x):
    """sin(x)"""
    return math.sin(x)


def main():
    print('\n', f1.__doc__)
    f1raiz1 = bisseccao(f1, intervalo=[1, 10])
    f1raiz2 = falsaposicao(f1, intervalo=[1, 10])
    f1raiz4 = pontofixo2(f1, x0=2.6)
    # This calculates the f1'
    f1_ = lambdify(x, diff(f1.__doc__))
    f1raiz6 = newtonraphson(f1, f1_, 2.6)

    f1raiz7 = secante(f1, 1, 10)

    print("Bisseccao", f1raiz1, f1(f1raiz1))
    print("FalsPosic", f1raiz2, f1(f1raiz2))
    print("PontoFix2", f1raiz4)
    print("NewtonRap", f1raiz6, f1(f1raiz6))
    print("Secante  ", f1raiz7, f1(f1raiz7))

    print('\n', f2.__doc__)
    f2raiz1 = bisseccao(f2, intervalo=[-1, 1])
    f2raiz2 = falsaposicao(f2, intervalo=[-1, 1])
    f2raiz3 = pontofixo1(f2, phi, x0=0.3)
    f2raiz4 = pontofixo2(f2, x0=0.3)
    f2raiz5 = optimize.fixed_point(f2, x0=0.3)

    # This calculates the f2'
    f2_ = lambdify(x, diff(f2.__doc__))

    f2raiz6 = newtonraphson(f2, f2_, 0.3)
    f2raiz7 = secante(f2, -1, 1)
    print("Bisseccao", f2raiz1, f2(f2raiz1))
    print("FalsPosic", f2raiz2, f2(f2raiz2))
    print("PontoFix1", f2raiz3, f2(f2raiz3))
    print("PontoFix2", f2raiz4, f2(f2raiz4))
    print("SCIPYMPF ", f2raiz5, f2(f2raiz5))
    print("NewtonRap", f2raiz6, f2(f2raiz6))
    print("Secante  ", f2raiz7, f2(f2raiz7))

    print('\n', f3.__doc__)
    f3raiz1 = bisseccao(f3, intervalo=[-1, 1], tol=eps)
    f3raiz2 = falsaposicao(f3, intervalo=[-1, 1], xtol=eps, ytol=eps)
    f3raiz4 = pontofixo2(f3, x0=0.1, tol=eps, maxiter=500)

    # This calculates the f3'
    f3_ = lambdify(x, diff(f3.__doc__))
    f3raiz6 = newtonraphson(f3, f3_, 0.5)
    f3raiz7 = secante(f3, -1, 1)

    print("Bisseccao", f3raiz1, f3(f3raiz1))
    print("FalsPosic", f3raiz2, f3(f3raiz2))
    print("PontoFix2", f3raiz4, f3(f3raiz4))
    print("NewtonRap", f3raiz6, f3(f3raiz6))
    print("Secante  ", f3raiz7, f3(f3raiz7))


if __name__ == "__main__":
    main()
