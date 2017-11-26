#!/usr/bin/env python3

import collections
import numbers
from numerico.geral.erros import eps


def falsaposicao(funcao, intervalo, xtol=eps, ytol=eps, maxiter=25):
    """Retorna a raiz de uma fucao em um intervalo dado

    Entrada:
    funcao: objeto: funcao a ser avaliada
    intervalo: iterable: intervalo fechado
    xtol: float: tolerancia do resultado de x
    ytol: float: tolerancia do resultado de f(x)
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
        x = ((a * f(b)) - b * f(a)) / (f(b) - f(a))

        if abs(b - a) < xtol:
            return x
        if abs(f(x)) < ytol:
            return x

        if f(x) * f(b) < 0:
            a = x
        else:
            b = x

    return x
