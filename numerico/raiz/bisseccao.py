#!/usr/bin/env python3

import collections
import numbers
from numerico.geral.erros import eps


def bisseccao(funcao, intervalo, tol=eps, maxiter=25):
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

    if not isinstance(intervalo, collections.abc.Sequence):
        raise TypeError("Intervalo deve ser iteravel")
    elif len(intervalo) != 2:
        raise ValueError("Intervalo deve possuir dois elementos")

    for n in intervalo:
        if not isinstance(n, numbers.Number):
            raise TypeError("Elementos do intervalo devem ser numericos")
    a = intervalo[0]
    b = intervalo[1]

    if f(a) * f(b) >= 0:
        msg = "Sinais da imagem devem ser diferentes f({})={} e f({}) = {}"
        raise ValueError(msg.format(a, f(a), b, f(b)))
    # END: Parameter checking }}}

    for iter in range(maxiter):
        x = (a + b) / 2

        if f(x) == 0:
            return x
        if abs(b - a) < tol:
            return x

        if f(x) * f(b) < 0:
            a = x
        else:

            b = x
    return x
