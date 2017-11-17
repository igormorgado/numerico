#!/usr/bin/env python3

import collections
import numbers

eps = 1e-7

def pontofixo1(funcao, phi, x0, xtol=eps, ytol=eps, maxiter=200):
    """Encontra o ponto onde funcao(x) == x

    Utiliza o algoritimo definido em Ruggiero pg. 64, onde a funcao phi
    do ponto fixo e' dada

    Entrada:
    funcao: objeto: funcao a ser avaliada
    x0: number: Ponto inicial
    tol: float: tolerancia do resultato
    maxiter: int: numero maximo de iteracoes

    Retorna:
    raiz: float: ponto no intervalo da raiz aproximada
    """

    # Parameter checking {{{
    if not callable(funcao):
        raise TypeError("Funcao deve ser executavel")
    f = funcao

    if not callable(phi):
        raise TypeError("phi deve ser executavel")

    if not isinstance(x0, numbers.Number):
        raise TypeError("Elementos do intervalo devem ser numericos")
    # END: Parameter checking }}}

    if abs(f(x0)) < ytol:
        return x0

    for iter in range(maxiter):
        x = phi(x0)

        if abs(f(x)) < ytol:
            return x

        if abs(x - x0) < xtol:
            return x

        x0 = x

    return x


def pontofixo2(funcao, x0, tol=eps, maxiter=200):
    """Encontra o ponto onde funcao(x) == x

    Entrada:
    funcao: objeto: funcao a ser avaliada
    x0: number: Ponto inicial
    tol: float: tolerancia do resultato
    maxiter: int: numero maximo de iteracoes

    Retorna:
    raiz: float: ponto no intervalo da raiz aproximada

    Usa o metodo de Steffensen baseando-se na aceleracao de convergencia de
    Aitken's. Veja: sBurden, Faires, "Numerical Analysis", 5 edicao, pg. 80
    """

    # Parameter checking {{{
    if not callable(funcao):
        raise TypeError("Funcao deve ser executavel")
    f = funcao

    if not isinstance(x0, numbers.Number):
        raise TypeError("Elementos do intervalo devem ser numericos")
    # END: Parameter checking }}}

    p0 = x0

    for iter in range(maxiter):
        p1 = f(p0)
        p2 = f(p1)
        d = p2 - 2 * p1 + p0

        if d == 0:
            return p2
        else:
            p = p0 - (p1 - p0) * (p1 - p0) / d

        if p0 == 0:
            relerr = p
        else:
            relerr = (p - p0) / p0

        if relerr < tol:
            return p

        p0 = p

    raise RunTimeError("Nao convergiu apos {} iteracoes {}.".format(iter, p))
