#!/usr/bin/python3

from sympy import symbols, Poly, Function, Symbol, expand, diff
from sympy.abc import x

x0, y0 = symbols('x0, y0')
x1, y1 = symbols('x1, y1')
x2, y2 = symbols('x2, y2')
x3, y3 = symbols('x3, y3')
x4, y4 = symbols('x4, y4')


# The lagrange coefficient
def lagrange_coefficient(x, i, X):

    l = 1
    for j in range(len(X)):
        if i != j:
            term = (x - X[j]) / (X[i] - X[j])
            l  =  l * term

    return l

def Lagrange_polynomial(x, X, Y):

    L = 0
    for i in range(len(X)):
        L += lagrange_coefficient(x, i, X) * Y[i] 

    return L


# Let the points
X = [x0, x1, x2, x3, x4]
Y = [y0, y1, y2, y3, y4]
L = Lagrange_polynomial(x, X, Y)
l = [ L.as_terms()[0][i][0] for i in range(len(X)) ]
