#!/usr/bin/env python3

import numbers
from numerico.geral.erros import import eps


def secante(funcao, x0, x1, xtol=eps, ytol=eps, maxiter=200):
    """Retorna a raiz de uma funcao utilizando o metodo das secantes

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

    if not isinstance(x0, numbers.Number):
        raise TypeError("Limite inferior do intervalo deve ser numerico")

    if not isinstance(x1, numbers.Number):
        raise TypeError("Limite superior do intervalo deve ser numerico")
    # END: Parameter checking }}}

    if abs(f(x0)) < ytol:
        return x0

    if abs(f(x1)) < ytol:
        return x1

    if abs(x1 - x0) < xtol:
        return x1

    for iter in range(maxiter):
        x2 = x1 - (f(x1) / (f(x1) - f(x0))) * (x1 - x0)

        if abs(f(x2)) < ytol:
            return x2

        if abs(x2 - x1) < xtol:
            return x2

        x0, x1 = x1, x2

    msg = "Funcao nao converge apos {} iteracoes: x={}"
    raise RuntimeError(msg.format(iter, x2))
