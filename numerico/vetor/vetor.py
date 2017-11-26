#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import randint
from random import random
from math import e as euler
from math import hypot
from math import cos, sin, acos, asin, sqrt, pi, atan
from numerico.geral.erros import eps, comozero
import types


# Quantas casas significativas mostrar em representacoes de string
casas = 4

def zeros_como(obj):
    """Retorna um objeto zerado do tipo passado"""

    if not (isinstance(obj, Vetor) or isinstance(obj, Matriz)):
        return NotImplemented

    return obj.zeros(*obj.dim())


def uns_como(obj):
    """Retorna um objeto unitario do tipo passado"""

    if not (isinstance(obj, Vetor) or isinstance(obj, Matriz)):
        return NotImplemented

    return obj.uns(*obj.dim())


class Vetor(object):
    """Classe para criar e manipular vetores"""

    def __init__(self, inputlist):
        """Initialize the new object, and verify initial conditions"""

        # Verify if the inputlist is a iterable but the first element isn't
        # it means is a single dimensional iterable
        try:
            self.n = len(inputlist)
        except TypeError:
            raise(VetorError, "Parametro de entrada nao possui tamanho")

        if hasattr(inputlist[0], "__iter__"):
            raise(VetorError, "Elementos nao sao escalares")

        self.elem = inputlist

    def __getitem__(self, index):
        """ Returns a element of a vector as scalar or a slice as a vector"""
        if isinstance(index, slice):
            return Vetor(self.elem[index])
        else:
            return self.elem[index]

    def __setitem__(self, index, value):
        """Set a vector element based on index"""
        self.elem[index] = value

    def __str__(self):
        """Returns a Printable Vetor Object"""
        return '[' + ', '.join([str(x) for x in self.elem]) + ']'

    def __repr__(self):
        """Returns the valid python representation of the Vetor"""
        modulename = str(type(self).__module__)

        ichars = len(str(int(self.max())))
        slen = ichars + casas
        fstr = "{{:>{}.{}g}}".format(slen, casas)

        if modulename == "__main__":
            s = str(type(self).__name__)
        else:
            s = modulename + '.' + str(type(self).__name__)

        s += '(['
        s += ', '.join([fstr.format(x) for x in self.elem])
        s += '])'

        return s

    def __len__(self):
        """Return the vector size (procedural)"""
        return self.n

    def len(self):
        """Return the vector size (object)"""
        return self.n

    def __add__(self, other):
        """ Sum two vectors and return a copy with result"""
        n = len(self)

        if n != len(other):
            raise(VetorError, "Vetor dimensions are not equal")

        v = zeros_como(self)

        for i in range(n):
            v[i] = self[i] + other[i]

        return v

    def __sub__(self, other):
        """Subtract two vectors and return a copy with result"""
        n = len(self)

        if n != len(other):
            raise(VetorError, "Vetor dimensions are not equal")

        v = zeros_como(self)

        for i in range(n):
            v[i] = self[i] - other[i]

        return v

    def __mul__(self, other):
        """Multiply self by other (if is a scalar) and return a copy of it.
        If other is also a vector returns element wise multiplication"""

        s = len(self)
        v = zeros_como(self)

        if isinstance(other, Vetor):
            # Both operands are Vetors
            # In this case perform a element wise product
            r = len(other)

            if s != r:
                raise(VetorError, "Vetor dimensions are not equal")

            for i in range(s):
                v[i] = self[i] * other[i]
        else:
            # check if other is a scalar
            if hasattr(other, "__len__"):
                raise(VetorError, "Operand isn't an scalar")

            for i in range(s):
                v[i] = other * self[i]

        return v

    def __rmul__(self, other):
        """Multiply vector by a scalar and return a copy of it."""

        s = len(self)
        v = zeros_como(self)

        # check if other is a scalar
        if hasattr(other, "__len__"):
            raise(VetorError, "Operand isn't an scalar")

        for i in range(s):
            v[i] = self[i] * other

        return v

    def __div__(self, other):
        """Divides self by a scalar and return a copy of it.
        If other is also a vector returns element wise division"""

        s = len(self)
        v = zeros_como(self)

        if isinstance(other, Vetor):
            # Both operands are Vetors
            # In this case perform a element wise product
            r = len(other)

            if s != r:
                raise(VetorError, "Vetor dimensions are not equal")

            for i in range(slen):
                v[i] = self[i] / float(other[i])
        else:
            # check if other is a scalar
            if hasattr(other, "__len__"):
                raise(VetorError, "Operand isn't an scalar")

            for i in range(s):
                v[i] = self[i] / float(other)

        return v

    def __rdiv__(self, scalar):
        """Divide a scalar by a vector is an Invalid operation"""
        raise(VetorError, "Not possible divide a scalar by a vector")

    def __eq__(self, other):
        """Returns True if two vectors are numericaly equal"""
        s = len(self)
        r = len(other)

        if s != r:
            raise(VetorError, "Vetor dimensions are not equal")

        # Two vectors are numericaly the same if the difference
        # between both of them are smaller than given precisao
        for i in range(s):
            if not comozero(self[i] - other[i]):
                return False

        return True

    def __abs__(self):
        """Returns vector absolute value"""
        v = zeros_como(self)

        for i in range(self.n):
            v[i] = abs(self[i])

        return v

    def __neg__(self):
        """Returns vector negative value"""
        v = zeros_como(self)

        for i in range(self.n):
            v[i] = -self[i]

        return v

    def __pow__(self, other):
        """Retorna a potencia de um vetor elemento a elemento"""
        n = len(self)

        v = zeros_como(self)

        for i in range(n):
            v[i] = self[i]**other

        return v

    @classmethod
    def _createVetor(cls, elem):
        """Metodo para criar uma instancia de vetor"""
        return cls(elem)

    @classmethod
    def zeros(cls, n):
        """Retorna um vetor nulo n-dimensional"""
        return cls._createVetor([0] * n)

    @classmethod
    def uns(cls, n):
        """Retorna um vetor unitario n-dimensional"""
        return cls._createVetor([1] * n)

    def index(self, value, start=0, stop=-1):
        """Retorna a primeira posicao de um elemento no vetor"""
        return self.elem.index(value, start, stop)

    def imin(self):
        """Retorna a posicao do menor elemento"""
        return self.index(min(self))

    def imax(self):
        """Retorna a posicao do maior elemento"""
        return self.elem.index(max(self))

    def min(self):
        """Retorna o menor elemento"""
        return min(self)

    def max(self):
        """Retorna o maior elemento"""
        return max(self)

    def dim(self):
        """Retorna a dimensao do vetor como uma u-pla"""
        return (self.n, )

    def pi(self, other):
        """Retorna o produto interno entre dois vetores"""

        s = len(self)
        r = len(other)

        if s != r:
            raise(VetorError, "Vetor dimensions are not equal")

        d = 0
        for i in range(s):
            d += self[i] * other[i]

        return d

    def pd(self, other):
        """Retorna o product diatico"""
        return Matriz([self]).T() * Matriz([other])

    def pv(self, other):
        """Retorna o produto vetorial de um vetor em R3"""

        assert self.n == other.n == 3, "Produto vetorial definido somente em R3"

        u, v = self, other

        return Vetor([u[1] * v[2] - u[2] * v[1],
                      u[2] * v[0] - u[0] * v[2],
                      u[0] * v[1] - u[1] * v[0]])

    def norma(self):
        """Retorna a norma euclidiana"""
        return (self.pi(self))**(0.5)

    def normaliza(self):
        """Retorna o vetor unitario"""
        return self * (1 / self.norma())

    def soma(self):
        """Retorna a soma dos elementos do vetor"""
        r = 0
        for i in range(len(self)):
            r += self[i]

        return r

    def prod(self):
        """Retorna o produto dos elementos do vetor"""
        r = 0
        for i in range(len(self)):
            r *= self[i]

        return r

    def conjugate(self):
        """Retorna o vetor conjugado"""
        v = zeros_como(self)
        for x in range(self.n):
            v[x] = (self[x]).conjugate()

        return v

    def e_ortogonal(self, other):
        """Verifica se o vetor e' ortogonal a outro"""
        if self.pi(other) == 0:
            return True
        else:
            return False

    def e_paralelo(self, other):
        """Verifica se o vetor e' paralelo a outro"""
        if (self == other) or (self.normaliza() == other.normaliza()):
            return True
        else:
            return False

    def polyeval(self, x):
        """Retorna o valor do vetor polinomial em 'x'
        Vetor([1,2,3]).polyeval(5) calcula 1x^2 + 2x + 3 com x=5
        """
        return NotImplemented

    def zeros_como(self):
        """Retorna um objeto zerado do tipo passado"""

        return self.zeros(*self.dim())

    def uns_como(self):
        """Retorna um objeto zerado do tipo passado"""

        return self.uns(*self.dim())

    def faca_zero(self):
        """Converte todos os valores menores que precisao para 0"""
        n = self.n

        for i in range(n):
            if comozero(self[i]):
                self[i] = 0

        return self

    def innulo(self):
        """Retorna o indice do primeiro elemento nao nulo do vetor"""
        for i in range(self.n):
            if not comozero(self[i]):
                return i
        return None


"""Classes para tratamento de erros"""

class VetorError(Exception):
    """Standard Error class for Vetor object"""
    pass


