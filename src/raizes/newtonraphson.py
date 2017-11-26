#!/usr/bin/env python3

import collections
import numbers
from erros import eps


def newtonraphson(funcao, funcao_, x0, xtol=eps, ytol=eps, maxiter=200):
    """Retorna a raiz de uma fucao em um intervalo dado

    Entrada:
    funcao: objeto: funcao a ser avaliada
    intervalo: list: limites do intervalo fechado
    tol: float: tolerancia do resultato
    maxiter: int: numero maximo de iteracoes

    Retorna:
    raiz: float: ponto no intervalo da raiz aproximada
    """

    # Parameter checking {{{
    if not callable(funcao):
        raise TypeError("Funcao deve ser executavel")
    f = funcao

    if not callable(funcao_):
        raise TypeError("funcao_ deve ser executavel")
    f_ = funcao_

    if not isinstance(x0, numbers.Number):
        raise TypeError("Chute inicial deve ser numerico")
    # END: Parameter checking }}}

    if abs(f(x0)) < ytol:
        return x0

    for iter in range(maxiter):
        x1 = x0 - f(x0) / f_(x0)

        if abs(f(x1)) < ytol:
            return x1

        if abs(x1 - x0) < xtol:
            return x1

        x0 = x1

    return x1
