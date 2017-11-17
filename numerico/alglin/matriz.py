#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .vetor import *

""" linalg.py implements common Linear Algebra structures and operation
    as studied in regular numerical linear algebra classes.

    Whatever simplicity vs speed decisions should be made this library
    will pickup the first, also it will avoid all Pythonisms for sake
    of explicit demonstration of algorithms.

    This library is distributed under GNU/GPLv2 and can be used for any
    purposes.

    author: Igor Morgado <morgado.igor@gmail.com>

    TODO: Enlarge the documentation explaining the mathematical theorems
        behind them and using latex annotations

        - Translate everything to portuguese
		- Implement more specific errors
		- Separate erros from main file (for sake of readibility)
		- Validation methods (istirangular, istriangsup, ...
		- Implement l squares
		- Implement tests for all matrices types and operations
        - Implementar dot matricial como um produto de matrizes. :-D
        - Implementar decomposicao LDU
        - Implementar LU compacta
        - Implementar diagonalizacao dominante
        - Implementar (em outro lugar) visualizacao animada das eliminacoes
          # plt.matshow(M.linhas) faz o trabalho, desenvolver uma funcao 
          # simplificada que atua direto na classe Matriz e Vetor.
        - Implementar metodos inplace/not inplace.
        - Reduces the ammount of object copies (lists, matrices, etc...)
        - Reduces mutations...
        - [ [0]*n for _ in range(n) ]  para uma matrix n x n.
        - O matematica representa o sistema Ax=b como uma lista com 2 elementos
          o primeiro elemento e A e o segundo b. Uma ideia, mas prefiro 
          faze-los separados. A uma matriz e B outra matriz assim posso juntar
          se quiser com AxB = [ A, B ] (boa ideia).
        - Metodos para converter matrizes em vetores e vetores em matrizes.
        - matriz2vetor = lambda x: Vetor(x[0].elem+x[1].elem)




    Good ideas:

    # Determinante by mms
    http://code.activestate.com/recipes/189971-basic-linear-algebra-matrix/

    # Get slice
    # House, ratval, qr, inverse, hessenberg, eigenvals, eigenvects
    http://users.rcn.com/python/download/python.htm

    # Funcao para criar matrizes triangulares superiores, inferiores e
    # diagonais aleatorias

    # Criar classes especificas para matrizes triangulares (superior e
    # inferior), diagonais, quadradas, tridiagonais, pentadiagonais 
    # de forma que seja possivel otimizar algumas operacoes e garantir algumas
    # propriedades algebricas.

    # Funcoes para criar matrizes de Hilbert

    # Uniformizar nomes de funcoes, ou usamos e_escaler ou eescalar ou eEscalar
    # mas nao os 3 em lugares diferentes.


"""

from __future__ import division
from random import randint
from random import random
from math import e as euler
from math import hypot
from math import cos, sin, acos, asin, sqrt, pi, atan
import types


# Quantas casas significativas mostrar em representacoes de string
casas = 4

# Valores abaixo de 'precisao' sao considerados como zero
precisao = 1e-9


def comozero(z, precisao=precisao):
    """Verifica se um valor e menor que a precisao numerica da maquina"""
    return abs(z) < precisao


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


def vandermode(v):
    """Retorna a matriz de Vandermode gerada por um vetor"""
    n = v.n
    M = Matriz.zeros(n)
    for j in range(n):
        M[j] = v**(n - j - 1)

    return M.T()


def e_escalar(n):
    """Verifica se o item e' um escalar"""
    return isinstance(n, types.IntType) or isinstance(n, types.FloatType) or isinstance(n, types.ComplexType)


def Rn(Ul):
    """ Retorna o produto de uma lista de transformacoes para formar R A^-1
    Se o retorno desta funcao eh multiplica do A a direita, o resultado e' R"""

    R = Matriz.identidade(Ul[0].m)

    for u in Ul:
        R = u * R

    return R


def Qn(Ul):
    """ Retorna o produto de uma lista de transformacoes para formar Q"""

    Q = Matriz.identidade(Ul[0].m)
    Ul.reverse()

    for u in Ul:
        Q = Q * u

    return Q.T()


class Matriz(object):
    """Classe para criar e manipular matrizes"""

    def __init__(self, inputlist):
        """Initializes the matrix from a iterable and verifies if it is a valid
        matrix"""

        if isinstance(inputlist, Vetor):
            self.m = 1
            self.n = len(inputlist)
            self.linhas = [inputlist]
        else:

            # Verify is the matrix is a bidimensional object
            # (and not anything else)
            try:
                self.m = len(inputlist)
            except TypeError:
                raise(MatrizError, "Matriz nao possui linhas")

            try:
                self.n = len(inputlist[0])
            except TypeError:
                raise(MatrizError, "Matriz nao possui colunas")

            # Verify if all lines has the same dimension
            for i in range(1, self.m):
                assert self.n == len(
                    inputlist[i]), "Numero de colunas nao e' homogeneo na linha {}".format(i)

            # TODO:  Each line need to be a treated as a Vetor object
            self.linhas = inputlist

    def __getitem__(self, index):
        """Retorna um valor da matriz, este e' um metodo padrao"""

        # How to implement slice...
        if isinstance(index, slice):
            return NotImplemented
        else:
            return Vetor(self.linhas[index])

    def __setitem__(self, index, value):
        """Altera o valor de uma matriz, este e' um metodo padrao"""

        # How to implement slice...
        if isinstance(index, slice):
            pass
        else:
            self.linhas[index] = value

    def __str__(self):
        """Exibe a matrix na forma linha coluna"""
        lines = []
        
        ichars = len(str(int(self.max())))
        slen = ichars + casas + 1
        fstr = "{{:>{}.{}g}}".format(slen, casas)
        
        for row in self.linhas:
            line = []
            for element in row:
                line.append(fstr.format(element))

            line = ' '.join(line)
            lines.append(line)

        lines = '\n'.join(lines)

        m, n = self.dim()
        rep = 'Matriz rank (%d, %d)\n%s\n' % (m, n, lines)

        return rep

    def __repr__(self):
        """Returns the valid python representation of the Matriz"""

        # Full qualified class name
        s = ''
        modulename = str(type(self).__module__)
        if modulename == "__main__":
            s += str(type(self).__name__)
        else:
            s += modulename + '.' + str(type(self).__name__)

        s += "([\n"
        skip = len(s) - 1

        lines = []

        # Isto acha o maior numero
        ichars = len(str(int(self.max())))
        slen = ichars + casas + 1
        fstr = "{{:>{}.{}g}}".format(slen, casas)

        for row in self.linhas:
            line = []
            for element in row:
                line.append(fstr.format(element))
            line = ', '.join(line)
            lines.append(' ' * skip + '[' + line + ']')

        lines = ',\n'.join(lines)

        s += lines
        s += "\n" + ' ' * (skip - 2) + "])"

        return s

    def __len__(self):
        """Retorna o numero de elementos da matriz"""
        return self.m * self.n

    def len(self):
        """Retorna o numero de elementos da matriz"""
        return self.m * self.n

    def __add__(self, other):
        """Adiciona a matriz atual a uma outra, retorna uma nova matriz"""

        if self.dim() != other.dim():
            raise(MatrizError, "Matrizes com dimensoes nao compativeis")

        # Define a matriz resultado
        r = zeros_como(self)

        for i in range(self.m):
            for j in range(self.n):
                r[i][j] = self[i][j] + other[i][j]

        return r

    def __sub__(self, other):
        """Subtracts two matrices and return a copy"""

        if self.dim() != other.dim():
            raise(MatrizError, "Matrizes com dimensoes nao compativeis")

        # Define a matriz resultado
        r = zeros_como(self)

        for i in range(self.m):
            for j in range(self.n):
                r[i][j] = self[i][j] - other[i][j]

        return r

    def __mul__(self, other):
        """Multiply self by other (handle scalar and matrix operations) and return a copy"""

        m, n = self.m, self.n

        if isinstance(other, type(self)):
            # Aqui estamos multiplicando duas matrizes
            om, on = other.m, other.n

            assert n == om, "Erro de erro de linhas e colunas entre matrizes"

            M = type(self).zeros(m, on)

            for i in range(m):
                for j in range(on):
                    M[i][j] = self.linha(i).pi(other.coluna(j))

        else:
            # Aqui esperamos multiplicar matriz por escalar
            M = zeros_como(self)

            for i in range(m):
                for j in range(n):
                    M[i][j] = other * self[i][j]

        return M

    def __rmul__(self, other):
        """Multiply self by other (handle scalar and matrix operations) and return a copy"""
        m, n = self.m, self.n
        M = zeros_como(self)

        for i in range(m):
            for j in range(n):
                M[i][j] = other * self[i][j]

        return M

    def __div__(self, other):
        """Divisao de matrizes nao e' bem definida"""
        return NotImplemented

    def __rdiv__(self, other):
        """Divisao de matrizes nao e' bem definida"""
        return NotImplemented

    def __eq__(self, other):
        """ Verifica se duas matrizes sao iguais numericamente """

        if self.dim() != other.dim():
            return False

        for i in range(self.m):
            for j in range(self.n):
                if not comozero(self[i][j] - other[i][j]):
                    return False

        return True

    def __abs__(self):
        """Returns the absolute value of matrix"""
        M = zeros_como(self)

        for i in range(self.m):
            for j in range(self.n):
                M[i][j] = abs(self[i][j])

        return M

    def __neg__(self):
        """Returns the negative value of matrix"""
        M = zeros_como(self)

        for i in range(self.m):
            for j in range(self.n):
                M[i][j] = -self[i][j]

        return M

    @classmethod
    def _criaMatriz(cls, inputlist):
        """Metodo interno para criar uma matriz """

        return cls(inputlist)

    @classmethod
    def zeros(cls, m=1, n=None):
        """Creates a zero matrix rank n,m"""

        if n is None: n = m

        assert m > 0, "Numero de linhas deve ser maior que zero"
        assert n > 0, "Numero de colunas deve ser maior que zero"

        M = []
        for i in range(m):
            M.append([0] * n)

        return cls._criaMatriz(M)

    @classmethod
    def uns(cls, m=1, n=None):
        """Creates a uns matrix rank n,m"""

        if n is None: n = m

        assert m > 0, "Numero de linhas deve ser maior que zero"
        assert n > 0, "Numero de colunas deve ser maior que zero"

        M = []
        for i in range(m):
            M.append([1] * n)

        return cls._criaMatriz(M)

    @classmethod
    def random(cls, m=1, n=None):
        """Cria uma matriz aleatoria de dimensoes m,n"""

        if n is None: n = m

        assert m > 0, "Numero de linhas deve ser maior que zero"
        assert n > 0, "Numero de colunas deve ser maior que zero"

        M = cls.zeros(m, n)
        for i in range(m):
            for j in range(n):
                M[i][j] = random()

        return M

    @classmethod
    def irandom(cls, m=1, n=None, interval=(-10, 10)):
        """Cria uma matriz aleatoria de inteiros com dimensoes m,n"""
        
        if n is None: n = m

        assert m > 0, "Numero de linhas deve ser maior que zero"
        assert n > 0, "Numero de colunas deve ser maior que zero"

        M = cls.zeros(m, n)
        for i in range(m):
            for j in range(n):
                M[i][j] = randint(interval[0], interval[1])

        return M

    @classmethod
    def identidade(cls, m=1):
        """Creates identidade matrix order n"""

        I = cls.zeros(m)

        assert m > 0, "Numero de linhas deve ser maior que zero"

        for i in range(m):
            I[i][i] = 1

        return I

    def max(self):
        """Retorna o maior elemento da matriz"""

        if len(self) == 0:
            return False
        else:
            m = self[0][0]

        for i in range(self.m):
            for j in range(self.n):
                if self[i][j] > m:
                    m = self[i][j]

        return m

    def min(self):
        """Retorna o menor elemento da matriz"""

        if len(self) == 0:
            return False
        else:
            m = self[0][0]

        for i in range(self.m):
            for j in range(self.n):
                if self[i][j] < m:
                    m = self[i][j]

        return m

    def conjugate(self):
        """Retorna a matriz conjugada"""
        M = zeros_como(self)
        for i in range(self.m):
            for j in range(self.n):
                M[i][j] = (self[i][j]).conjugate()

        return M

    def dim(self):
        """Retorna as dimensoes da matriz"""
        return self.m, self.n

    def coluna(self, idx):
        """Retorna a coluna definida no indice"""
        return Vetor([l[idx] for l in self.linhas])

    def linha(self, idx):
        """Retorna a coluna definida no indice"""
        return Vetor(self[idx])

    def det(self):
        """Calcula o determinante da matriz"""
        M = self

        assert M.m == M.n, "A matriz deve ser quadrada"

        # Definicao de determinante 2x2
        if M.m == 2:
            return M[0][0] * M[1][1] - M[0][1] * M[1][0]

        return M.cofator_linha(0)

    def cofator_linha(self, m):
        """Calcula o determinante pela expansao de cofatores """

        M = self
        assert (m < M.m)
        d = 0
        for j in range(M.n):
            d += (-1)**(m+j) *  M[m][j] * M.submatriz(m,j).det() 

        return d

    def submatriz(self, m=0, n=0):
        """Retorna uma submatriz sem a linha m e coluna n"""
        M = type(self).zeros(self.m - 1, self.n - 1)

        # x: contador de linhas de M
        # y: contador de colunas de M
        x = 0

        for i in range(self.m):
            y = 0
            if (i == m):
                continue
            for j in range(self.n):
                if (j == n):
                    continue
                M[x][y] = self[i][j]
                y += 1
            x += 1

        return M

    def vdiv(self, pos):
        """Separa verticalmente a matriz, retornando uma lista com duas matrizes
        Entrando a matriz M retorna duas matrizes A, B
                                        [ A ]
                                        -----
                                        [ B ]
        """

        m, n = self.m, self.n
        M = self

        assert pos is not None, "Deve ser passada uma posicao para divisao"
        assert pos > 0, "Posicao de corte deve ser a partir da primeira posicao" 
        assert pos < m, "Posicao de corte deve ser antes do fim da matriz"

        Am, An =   pos, n
        Bm, Bn = m-pos, n

        A = self.zeros(Am, An)
        B = self.zeros(Bm, Bn)

        for i in range(m):
            for j in range(n):
                if i < pos:
                    A[i][j] = M[i][j]
                else:
                    B[i-pos][j] = M[i][j]


        return A, B

    def hdiv(self, pos):
        """Separa horizontalmente a matriz, retornando uma lista duas matrizes
        Entrando a matriz M retorna duas matrizes  [ A | B ]
        """

        m, n = self.m, self.n
        M = self

        assert pos is not None, "Deve ser passada uma posicao para divisao"
        assert pos > 0, "Posicao de corte deve ser a partir da primeira posicao" 
        assert pos < n, "Posicao de corte deve ser antes do fim da matriz"

        Am, An = m, pos
        Bm, Bn = m, n-pos

        A = self.zeros(Am, An)
        B = self.zeros(Bm, Bn)

        for i in range(m):
            for j in range(n):
                if j < pos:
                    A[i][j] = M[i][j]
                else:
                    B[i][j-pos] = M[i][j]

        return A, B

    def posto(self):
        """Retorna o posto de uma matriz. Isto e o numero de vetores linha
        nao nulo"""

        M = self.copia()
        M = M.rref()

        p=0
        for i in range(M.m):
            if M[i].innulo() is not None:
                p += 1

        return p
    
    def norma(self, tipo='inf'):
        """Retorna a norma matricial
        Por padrao retorna a norma do infinito, mas pode retornar a norma
        euclidiana ou um"""

        def norma_inf(m):
            s = []
            for i in range(self.m):
                s.append(abs(self.linha(i)).soma())

            return max(s)

        def norma_um(m):
            s = []
            for j in range(self.n):
                s.append(abs(self.coluna(j)).soma())

            return max(s)

        def norma_dois(m):
            s = 0
            for i in range(self.m):
                for j in range(self.n):
                    s += self[i][j]**2

            return s**(1 / 2.)

        my_norma = {
            'inf': norma_inf,
                'um': norma_um,
                'dois': norma_dois
        }

        try:
            return my_norma[tipo](self)
        except:
            return NotImplemented

    def cond(self):
        """Analisa o condicionamento de uma matriz
        cond(A) = ||A|| * ||A^-1||
        Se cond(A) -> 1 entao a matriz e bem condicionada.
        """
        return NotImplemented

    def zeros_como(self):
        """Retorna um objeto zerado do tipo passado"""
        return self.zeros(*self.dim())

    def uns_como(self):
        """Retorna um objeto unitario do tipo passado"""
        return self.uns(*self.dim())

    def aumenta(self, other):
        """ Retorna a matriz A aumentada da B como em [A | B ].
        O numero de linhas em ambas matrizes deve ser o mesmo. """

        if isinstance(other, Vetor):
            other = type(self)(other).T()

        m, n = self.m, self.n
        om, on = other.m, other.n

        assert m == om, "O numero de linhas das matrizes deve ser o mesmo"

        M = type(self).zeros(m, n + on)
        Mm, Mn = M.m, M.n

        for i in range(Mm):
            for j in range(Mn):
                if j < n:
                    M[i][j] = self[i][j]
                else:
                    # Pegue da matriz B
                    M[i][j] = other[i][j - n]

        return M

    def concat(self, other):
        """ Retorna a matriz A concatenada da B como em
           [ A ]
           [ B ]

        O numero de colunas em ambas matrizes deve ser o mesmo."""

        if isinstance(other, Vetor):
            other = type(self)([other])

        m, n = self.m, self.n
        om, on = other.m, other.n

        assert n == on, "O numero de colunas das matrizes deve ser o mesmo"

        M = type(self).zeros(m + om, n)
        Mm, Mn = M.m, M.n

        for i in range(Mm):
            for j in range(Mn):
                if i < m:
                    M[i][j] = self[i][j]
                else:
                    # Pegue da matriz B
                    M[i][j] = other[i - m][j]

        return M

    def T(self):
        """Retorna a transposta da matriz"""
        M = type(self).zeros(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                M[i][j] = self[j][i]

        return M

    def H(self):
        """Retorna a matriz adjunta ou conjugada transposta"""
        M = self.T()
        return M.conjugate()

    def e_quadrada(self):
        return self.m == self.n

    def diag(self):
        """Matriz Diagonal"""
        assert self.e_quadrada(), "Matriz nao e' quadrada"

        v = []
        for x in range(self.m):
            v.append(self[x][x])

        return Vetor(v)

    def traco(self):
        """Retorna o traco da matriz"""
        return self.diag().soma()

    def copia(self):
        """Retorna uma copia da matriz"""

        M = zeros_como(self)

        for i in range(self.m):
            for j in range(self.n):
                M[i][j] = self[i][j]

        return M

    def xgauss(self, inferior=False):
        """Gerador para o algoritimo de Gauss
        Se inferior=True tambem retorna a matriz dos escalares"""

        m, n = self.m, self.n

        M = self.copia()
        
        # Se solicitada a matriz inferior de escalares
        # cria a area para tal
        if inferior: 
            C = type(self).identidade(m)

        # Criando a matriz triangular superior
        for k in range(n):
            for i in range(k + 1, m):
                if not comozero(M[i][k]):
                    # Calcula lambda
                    lam = M[i][k] / M[k][k]

                    if inferior:
                        R = M.copia(), C.copia()
                        C[i][k] = lam
                    else:
                        R = M.copia()
                    yield R
                    M[i] = M[i] - lam * M[k]

        if inferior:
            R = M.copia(), C.copia()
        else:
            R = M.copia()

        yield R

    def gauss(self,inferior=False):
        """Retorna o ultimo elemento do gerador de Gauss"""
        for M in self.xgauss(inferior): pass
        
        return M

    def pivorev(self, m=None, n=None):
        """Retorna a posicao do pivo reverso comecando de m,n (indexado em 0)
        Se m,n nao for especificado usa a ultima posicao como inicio"""

        if m is None or m > self.m-1: m = self.m-1
        if n is None or n > self.n-1: n = self.n-1

        M = self

        for i in range(m,-1,-1):
            for j in range(n+1):
                if not comozero(M[i][j]):
                    return i,j

    def xjordan(self):
        """Gerador para o algoritimo de Jordan"""
        # BUG: Para matrizes LD. xrref, funciona.
        while True:
            yield "O uso de xjordan() esta desabilitado no momento"
            
        m, n = self.m, self.n

        M = self.copia()

        pi, pj = M.pivorev()

        # Criando a matriz triangular inferior
        # Se aplicado o algoritimo de Gauss antes
        # retorna uma matriz diagonal

        # Inicia da k-esima linha onde foi encontrado o pivo inferior.
        for k in range(pi,-1,-1):
            # Verifica se a diagonal da linha e zero.
            if comozero(M[k][k]):
                continue

            # Torna o pivo da linha atual unitaria
            M[k] = M[k] * (1/M[k][k])
            # Avalia todas as linhas acima da atual
            for i in range(k-1,-1,-1):
                yield M.copia()
                if not comozero(M[i][k]):
                     lam = M[i][k] / M[k][k]
                     M[i] = M[i] - lam * M[k]
                
        yield M.copia()
    
    def jordan(self):
        """Retorna o ultimo elemento do gerador de Jordan"""
        for M in self.xjordan(): pass
        
        return M

    def xgaussjordan(self):
        """Retorna todos os passos da eliminacao de Gauss-Jordan"""

        for g in self.xgauss(): yield g
        for j in g.xrref(): yield j

    def gaussjordan(self):
        """Executa a eliminacao de Gauss-Jordan na matriz dada"""
        return self.gauss().rref().faca_zero()

    def xrref(self):
        """Gerador para o algoritimo de Jordan"""
        m, n = self.m, self.n

        M = self.copia()

        # Inicia da k-esima linha onde foi encontrado o pivo inferior.
        for k in range(m-1,-1,-1):
            p = M[k].innulo()

            # Linha vazia
            if p is None: continue

            # Torna o pivo da linha atual unitaria
            M[k] = M[k] * (1/M[k][p])
            # Avalia todas as linhas acima da atual
            for i in range(k-1,-1,-1):
                yield M.copia()
                if not comozero(M[i][p]):
                     lam = M[i][p] / M[k][p]
                     M[i] = M[i] - lam * M[k]
                
        yield M.copia()

    def rref(self):
        """Retorna o ultimo elemento do gerador de rref"""
        for M in self.xrref(): pass
        
        return M

    def substreversa(self):
        """Aplica a substituicao reversa de uma matriz [ A | b ]
        Retorna o a matriz e o vetor solucao""" 
        # Talvez seja melhor passar A e b de forma separada
        # ao inves de [ A | b ].
        # Isso pode resolver um monte de bugs inclusive em outras funcoes.

        m, n = self.m, self.n

        M = self

        # Vetor resposta (fazer funcionar com B sendo uma matriz)
        x = Vetor.zeros(m)

        # Por enquanto b e' sempre a ultima coluna
        # Pode ser que no futuro implemente B como matriz e
        # multiplas solucoes
        b = M.coluna(n - 1)

        # n-2 aqui pq b faz parte da matriz.
        # do contrario seria n-1.
        for k in range(n - 2, -1, -1):
            # Aqui e n-q pq b faz parte da pmatriz
            # do contrario seria n
            x[k] = (b[k] - M[k][k:n - 1].pi(x[k:n - 1])) / M[k][k]

        return x

    def substreversa2(self, b):
        """Aplica a substituicao reversa de uma matriz [ A | b ]
        Retorna o a matriz e o vetor solucao""" 
        # Talvez seja melhor passar A e b de forma separada
        # ao inves de [ A | b ].
        # Isso pode resolver um monte de bugs inclusive em outras funcoes.

        m, n = self.m, self.n

        M = self

        # Vetor resposta (fazer funcionar com B sendo uma matriz)
        x = Vetor.zeros(m)

        # Por enquanto b e' sempre a ultima coluna
        # Pode ser que no futuro implemente B como matriz e
        # multiplas solucoes
        b = M.coluna(n - 1)

        # n-2 aqui pq b faz parte da matriz.
        # do contrario seria n-1.
        for k in range(n - 2, -1, -1):
            # Aqui e n-q pq b faz parte da pmatriz
            # do contrario seria n
            x[k] = (b[k] - M[k][k:n - 1].pi(x[k:n - 1])) / M[k][k]

        return x

    def substprogressiva(self):
        """Aplica a substituicao progressiva em uma matriz do tipo [A|B] """
        m, n = self.m, self.n

        M = self

        # Vetor resposta (fazer funcionar com B sendo uma matriz)
        x = Vetor.zeros(m)

        # Por enquanto b e' sempre a ultima coluna
        # Pode ser que no futuro implemente B como matriz e
        # multiplas solucoes
        b = M.coluna(n - 1)

        # n-2 aqui pq b faz parte da matriz.
        # do contrario seria n-1.
        for i in range(n-1):
            # Aqui e n-q pq b faz parte da pmatriz
            # do contrario seria n
            x[i] = (b[i] - M[i][:i+1].pi(x[i+1])) / M[i][i]

        return x

    def elimgauss(self):
        """Retorna a solucao da aplicando o algoritimo de
        Gauss + Substituicao Reversa"""

        # Aplica a reducao de Gauss
        M = self.gauss()

        # Aplica a substituicao reversa
        x = M.substreversa()

        return M,x

    def pivoteamento(self):
        """Retorna a matriz de pivoteamento"""
        m, n = self.m, self.n

        p = list(range(m))

        # Matriz nula quadrada
        P = type(self).zeros(m)

        M = self.copia().tamanhorelativo()

        for k in range(m):
            x = M.coluna(k)[k:].imax() 

            # P e' uma posicao relativa entao soma ao contador
            # pois a posicao e' sempre consumida
            x += k 

            # Melhor candidato nao e' a linha atual
            if x != k:
                # Entao troca linha
                M[k], M[x] = M[x], M[k]

                # Troca posicao da lista de trocas
                p[k], p[x] = p[x], p[k]

                
        # Retorna a matriz de pivoteamento
        #P = zeros_como(self)

        for i in range(m): 
            P[i][p[i]] = 1
        
        return P

    def decompLU(self):
        """Realiza decomposicao LU para uma matriz dada retorna as matrizes:
        P L U , onde P e' uma matriz de pivoteamento, L triangular inferior e 
        U triangular superior"""

        A = self

        m, n = A.m, A.n
    
        P = A.pivoteamento()

        A = P * A

        U, L = A.gauss(inferior=True)

        return P, L.faca_zero(), U.faca_zero()

    def solLU(self, b, metodo=None):
        """Retorna a solucao utilizando o algoritimo LU com pivoteamento
        Seja o sistema de eq. Lineares
            Ax=b

        Podemos decompor a matriz A em duas outras matrizes L, U (com
        pivoteamento parcial). Assim temos:

            PA = LU
            
            LUx = P b^T

        Assim ireimos primeiro solucionar o sistema (por subst progressiva):

            Ly = Pb

        Depois vamos resolver por subst. regressiva o sistema
        
            Ux = y
        
        Metodos podem ser:
            Crout:  Lower-Upper (Crout)
            LDU: Lower-Diagonal-Upper
            PLU: Pivoting Lower Upper
            PLDU: Pivoting Lower Diagonal Upper
            Cholesky: Cholesky Lower-Upper
            QRGS: QR por gram-schmidt
            SVD: Decomposicao de valores singulares
        """

        if self.e_simetrica():
            L,U = self.cholesky()
        else:
            P,L,U = self.decompLU()
            L = P*L

        n = L.n
        y = Vetor.zeros(n)
        x = Vetor.zeros(n)

        # Passo 1, soluciona L y = b
        y[0] = b[0]/L[0][0]
        for i in range(1,n):
            s = 0.0
            for j in range(i):
                s += L[i][j] * y[j]

            y[i] = (b[i] - s)/L[i][i]


        # Passo 2, soluciona Ux = y (regressivo)
        x[n-1] = y[n-1]/U[n-1][n-1]
        for i in range(n-1,-1,-1):
            s = y[i]
            for j in range(i+1, n):
                s = s - U[i][j] * x[j]
            x[i] = s/U[i][i]

        return x

        
    def LDU(self):
        """Implementa a decomposicao LDU"""
        return NotImplemented

    def gs(self):
        """Aplica o processo de ortogonalizacao de Gram-Schimidt na matriz
        vetor linha dada""" 

        B=self
        A=zeros_como(B)

        for i in range(B.m):
            A[i] = B[i]
            for k in range(i):
                A[i] -= ( B[i].pi(A[k]) / A[k].pi(A[k]) ) * A[k]

        return A

    def cholesky(self):
        """Aplica a decomposicao de Cholesky retornando L"""

        M = self.copia()

        L = zeros_como(M)

        for k in range(L.m):
            for i in range(0,k):
                s = 0
                for j in range(0,i):
                    s += M[i][j] * M[k][j]

                M[k][i] = (M[k][i] - s)/M[i][i]

            s = 0
            for j in range(0,k):
                s += M[k][j] * M[k][j]
            M[k][k] = (M[k][k] - s)**(.5)

        for k in range(L.m):
            L[k][0:k+1] = M[k][0:k+1]

        return (L, L.T())

    def QRGS(self):
        """Implementa a decomposicao QR por Gram-Schmidt"""

        M = self.copia()
        Mgs = M.gs()        # Aplica GS
        Q = Mgs.normaliza()  # Orto normal (Q)
        Q[0] *= -1
        R = Q * M
        QT = Q.T()

        #  Q * R = M
        return QT, R

    def QRREF(self):
        """QR por reflexoes (de Household
        Baseado em http://math.nist.gov/tnt/"""
        from numpy import sqrt

        QR = self.copia()
        m, n = QR.m, QR.n
        mm = min(m,n)
        Rd = Vetor.zeros(n)
        Q = Matriz.zeros(m,n)
        R = Matriz.zeros(n,n)

        for l in range(mm):
            norm = 0.0
            for linha in range(l, m):
                norm += QR[linha][l] * QR[linha][l]
            a = sqrt(norm)
            if QR[l][l] > 0:
                a = -a
            Rd[l] = a
            #Rd[l] = -(abs(QR[l][l:].norma()))

            if a != 0.0:
                QR[l][l] -= Rd[l]
                for col in range(l+1,n):
                    p = 0.0
                    for linha in range(l, m):
                        p -= QR[linha][col] * QR[linha][l]
                    p /= Rd[l] * QR[l][l]
                    for linha in range(l, m):
                        QR[linha][col] -= p * QR[linha][l]

        # R
        for linha in range(mm-1, -1, -1):
            R[linha][linha] = Rd[linha]
            for col in range(linha + 1, n):
                R[linha][col] = QR[linha][col]

        # Q
        for l in range(m-1, mm -1, -1):
            Q[l][l] = 1.0
        for l in range(mm-1, -1, -1):
            Q[l][l] = 1.0
            if QR[l][l] != 0.0:
                for col in range(l, m):
                    p = 0.0
                    for linha in range(l,m):
                        p -= Q[linha][col] * QR[linha][l]
                    p /= Rd[l] * QR[l][l]

                    for linha in range(l, m):
                        Q[linha][col] -= p * QR[linha][l]

        return Q, R

    def G(self, p, q):
        """ Retorna a matriz de rotacao que zera o elemento M[p][q](atan)"""

        k = self.m

        assert(p < k), "Eixo de rotacao 'p' deve ser menor que 'dim(M)'"
        assert(q < k), "Eixo de rotacao 'q' deve ser menor que 'dim(M)'"
        assert(p != q), "Eixos de rotacao devem ser diferentes"
        
        g = type(self).identidade(k)

        if self[p][q] != 0: 
            teta = atan( -self[p][q] / self[q][q] )
            g[p][p] = g[q][q] = cos(teta)
            g[p][q] = sin(teta)
            g[q][p] = -g[p][q]

        return g

    def U(self, p, q):
        """ Retorna a matriz de rotacao que zera o elemento M[p][q](sqrt)"""

        k = self.m

        assert(p < k), "Eixo de rotacao 'p' deve ser menor que 'dim(M)'"
        assert(q < k), "Eixo de rotacao 'q' deve ser menor que 'dim(M)'"
        assert(p != q), "Eixos de rotacao devem ser diferentes"
        
        u = type(self).identidade(k)

        # Onde esta pq era qp onde esta qq era pp ?
        # Verifica se retorna o mesmo que G para todo G
        if self[p][q] != 0: 
            norma = sqrt(self[q][q]**2 + self[p][q]**2)
            c = self[q][q]/norma
            s = self[p][q]/norma
            u[p][p] = u[q][q] = c
            u[p][q] = s
            u[q][p] = -s

        return u

    def Glist(self):
        """ Retorna uma lista de transformacoes de U que zeram abaixo da diagonal
        Executa de baixo para cima, esq para direita """

        T = self.copia()
        Gl = []

        # Zera somente abaixo da diagonal
        # e cria a lista de transformacoes Gl
        for j in range(T.n):
            for i in range(T.m-1,j,-1):
                Gt = T.G(i,j)
                T = Gt * T
                Gl.append(Gt)

        return Gl

    def Ulist(self):
        """ Retorna uma lista de transformacoes de U que zeram abaixo da diagonal
        cima para baixo, esq para direita """

        T = self.copia()
        Ul=[]

        # Zera somente abaixo da diagonal
        for j in range(T.n):
            for i in range(j+1,T.m):
                Ut = T.U(i,j)
                T = Ut * T
                Ul.append(Ut)

        return Ul

    def QRROT(self):
        """Retorna a decomposicao QR de A por rotacoes"""

        Gl = self.Glist()
        R = Rn(Gl) * self
        Q = Qn(Gl)

        return Q, R

    def QRrotU(self):
        """Retorna a decomposicao QR de A por rotacoes(using Ulist)"""

        Gl = self.Ulist()
        R = Rn(Gl) * self
        Q = Qn(Gl)

        return Q, R

    def QRfrancis(self, tol=1e-10, iters=10):
        """Retorna autovalores  QR por rotacoes"""
        # VERIFICAR APOS CORRECAO DA FUNCAO Qn e Rn

        T = self.copia()

        for i in range(iters):
            Ul = self.Ulist()
            R = Rn(Ul) * T
            Q = Qn(Ul)
            T = R * Q

            if ( erroFrancis(T) < tol ):
                break

        return Q,R

    def erroFrancis(self):
        """ Retorna o valor absoluto do maior elemento abaixo da diagonal"""

        maior = abs(self[1][0])
        for j in range(self.n):
            for i in range(j+1,self.m):
                if abs(self[i][j]) > maior: maior = abs(self[i][j])

        return maior

    def tetapq(self, p, q):
        """Retorna teta que elimina o elemento M[p][q] para uma matriz de rotacao"""

        assert(self.m == self.n), "Matriz deve ser quadrada"
        assert(p < self.m), "Eixo de rotacao 'p' deve ser menor que 'k'"
        assert(q < self.m), "Eixo de rotacao 'q' deve ser menor que 'k'"
        assert(p != q), "Eixos de rotacao devem ser diferentes"
        
        return acos(self[p][p]/(sqrt(self[p][p]**2 + self[q][p]**2)))

    @classmethod
    def matrot(cls, teta, k=2, p=0, q=None):
            """Retorna a matriz de rotacao k-dimensional nos eixos p e q de angulo
            teta"""

            if q is None: q = p+1

            assert(p < k), "Eixo de rotacao 'p' deve ser menor que 'k'"
            assert(q < k), "Eixo de rotacao 'q' deve ser menor que 'k'"
            assert(p != q), "Eixos de rotacao devem ser diferentes"

            U = cls.identidade(k)

            U[p][p] = U[q][q] = cos(teta)
            U[q][p] = sin(teta)
            U[p][q] = -U[q][p]

            return U

    def solveQR(self, b, metodo='GS'):
        """ Resolve usando QR, Rx = Q^T b  (default GS) """

        M = self.copia()
        n = M.n
        x = Vetor.zeros(n)

        if metodo == 'GS': 
            Q, R = M.QRGS()
        elif metodo == 'ROT':
            Q, R = M.QRROT()
        elif metodo == 'REF':
            Q, R = M.QRREF()
        else:
            return NotImplemented

        # Passo 1, resolve y = Q^T * b
        #y = Q.T() * Matriz(b).T()
        y = Matriz(b) * Q

        # Passo 2, resolve o sistema em R (Rx = y)
        x[n-1] = y[0][n-1] * (1/R[n-1][n-1])
        x[n-2] = ( y[0][n-2] - R[n-2][n-1] * x[n-1]) * (1/R[n-2][n-2])
        for j in range(n-3, -1, -1):
            s = 0.0
            for k in range(j+1, n):
                s += R[j][k] * x[k]
            x[j] = ( y[0][j] - s ) * (1/R[j][j])

        return x

    def SVD(self):
        """Retornar a decomposicao S,V,D"""
        
        A = self.copia()
        m,n = A.m , A.n
        nu = min(m, n)

        ncol, nlin = min(m-1, n), max(0, min(n-2, m))
        S = [[0.0 for i in range(nu)] for j in range(m)]  # singular values
        V = [0.0 for i in range(min(m + 1, n))]
        D = [[0.0 for i in range(n)] for j in range(n)]
        e = [0.0 for i in range(n)]
        sigma = [0.0 for i in range(m)]

        # Diagonalize the original matrix in superdiagonal V
        for k in range(max(ncol, nlin)):
            km_range = range(k, m)
            kp1_range = range(k + 1, n)
            if k < ncol:
                # Compute the 2-norm of k-th column for S
                V[k] = 0.0
                for i in km_range:
                    V[k] = hypot(V[k], A[i][k])
                if V[k] != 0.0:
                    if A[k][k] < 0.0:
                        V[k] = -V[k]
                    for i in km_range:
                        A[i][k] = A[i][k] / V[k]
                    A[k][k] = A[k][k] + 1.0
                V[k] = -V[k]
            for j in kp1_range:
                if (k < ncol) and (V[k] != 0.0):
                    t = 0.0
                    for i in km_range:
                        t = t + A[i][k] * A[i][j]
                    t = -t / A[k][k]
                    for i in km_range:
                        A[i][j] = A[i][j] + t * A[i][k]
                    e[j] = A[k][j]
            if k < ncol:
                for i in km_range:
                    S[i][k] = A[i][k]

            if k < nlin:
                # Compute the 2-norm of k-th column for D
                e[k] = 0.0
                for i in kp1_range:
                    e[k] = hypot(e[k], e[i])
                if e[k] != 0.0:
                    if e[k + 1] < 0.0:
                        e[k] = -e[k]
                    for i in kp1_range:
                        e[i] = e[i] / e[k]
                    e[k + 1] = e[k + 1] + 1.0
                e[k] = -e[k]
                if (k + 1 < m) and (e[k] != 0.0):
                    for i in range(k + 1, m):
                        sigma[i] = 0.0
                    for j in kp1_range:
                        for i in range(k + 1, m):
                            sigma[i] = sigma[i] + e[j] * A[i][j]
                    for j in kp1_range:
                        t = -e[j] / e[k + 1]
                        for i in range(k + 1, m):
                            A[i][j] = A[i][j] + t * sigma[i]
                for i in kp1_range:
                    D[i][k] = e[i]

        # setup bidiagonal matrix
        p = min(n, m + 1)
        if ncol < n:
            V[ncol] = A[ncol][ncol]
        if m < p:
            V[p - 1] = 0.0
        if (nlin + 1) < p:
            e[nlin] = A[nlin][p - 1]
        e[p - 1] = 0.0

        # generate S
        for j in range(ncol, nu):
            for i in range(m):
                S[i][j] = 0.0
            S[j][j] = 1.0
        for k in range(ncol - 1, -1, -1):
            if V[k] != 0.0:
                for j in range(k + 1, nu):
                    t = 0.0
                    for i in range(k, m):
                        t = t + S[i][k] * S[i][j]
                    t = -t / S[k][k]
                    for i in range(k, m):
                        S[i][j] = S[i][j] + t * S[i][k]
                for i in range(k, m):
                    S[i][k] = -S[i][k]
                S[k][k] = 1.0 + S[k][k]
                for i in range(k - 1):
                    S[i][k] = 0.0
            else:
                for i in range(m):
                    S[i][k] = 0.0
                S[k][k] = 1.0

        # generate D
        for k in range(n - 1, -1, -1):
            if (k < nlin) and (e[k] != 0.0):
                for j in range(k + 1, nu):
                    t = 0.0
                    for i in range(k + 1, n):
                        t = t + D[i][k] * D[i][j]
                    t = -t / D[k + 1][k]
                    for i in range(k + 1, n):
                        D[i][j] = D[i][j] + t * D[i][k]
            for i in range(n):
                D[i][k] = 0.0
            D[k][k] = 1.0

        # SVD
        pp = p - 1
        itr = 0
        eps = euler ** (-64)
        while p > 0.0:
            for k in range(p - 2, -2, -1):
                if k == -1:
                    break
                if abs(e[k]) <= eps * (abs(V[k]) + abs(V[k + 1])):
                    e[k] = 0.0
                    break
            if k == (p - 2):
                caso = 4
            else:
                for ks in range(p - 1, k - 1, -1):
                    if ks == k:
                        break
                    if ks != p:
                        t1 = abs(e[ks])
                    else:
                        t1 = 0.0
                    if ks != (k + 1):
                        t2 = abs(e[ks - 1])
                    else:
                        t2 = 0.0
                    t = t1 + t2
                    if abs(V[ks]) <= eps * t:
                        V[ks] = 0.0
                        break
                if ks == k:
                    caso = 3
                elif ks == (p - 1):
                    caso = 1
                else:
                    caso = 2
                    k = ks
            k += 1

            if caso == 1:
                f = e[p - 2]
                e[p - 2] = 0.0
                for j in range(p - 2, k - 1, -1):
                    t = hypot(V[j], f)
                    cs = V[j] / t
                    sn = f / t
                    V[j] = t
                    if j != k:
                        f = -sn * e[j - 1]
                        e[j - 1] = cs * e[j - 1]
                    for i in range(n):
                        t = cs * D[i][j] + sn * D[i][p - 1]
                        D[i][p - 1] = -sn * D[i][j] + cs * D[i][p - 1]
                        D[i][j] = t
            elif caso == 2:
                f = e[k - 1]
                e[k - 1] = 0.0
                for j in range(k, p):
                    t = hypot(V[j], f)
                    cs = V[j] / t
                    sn = f / t
                    V[j] = t
                    f = -sn * e[j]
                    e[j] = cs * e[j]
                    for i in range(m):
                        t = cs * S[i][j] + sn * S[i][k - 1]
                        S[i][k - 1] = -sn * S[i][j] + cs * S[i][k - 1]
                        S[i][j] = t
            elif caso == 3:
                # qr-elimination
                s = max(abs(V[k]), abs(e[k]), abs(e[p - 2]))
                s = max(s, abs(V[p - 2]), abs(V[p - 1]))
                memv = (V[p - 1], V[p - 2], e[p - 2], V[k], e[k])
                (sp, spm1, epm1, sk, ek) = map((lambda x: x / s), memv)
                b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0
                c = (sp * epm1) * (sp * epm1)
                shift = 0.0
                if (b != 0.0) or (c != 0.0):
                    shift = hypot(b, c)
                    if b < 0.0:
                        shift = -shift
                    shift = c / (b + shift)
                f = (sk + sp) * (sk - sp) + shift
                g = sk * ek
                for j in range(k, p - 1):
                    t = hypot(f, g)
                    cs = f / t
                    sn = g / t
                    if j != k:
                        e[j - 1] = t
                    f = cs * V[j] + sn * e[j]
                    e[j] = cs * e[j] - sn * V[j]
                    g = sn * V[j + 1]
                    V[j + 1] = cs * V[j + 1]
                    for i in range(n):
                        t = cs * D[i][j] + sn * D[i][j + 1]
                        D[i][j + 1] = -sn * D[i][j] + cs * D[i][j + 1]
                        D[i][j] = t
                    t = hypot(f, g)
                    (cs, sn) = (f / t, g / t)
                    V[j] = t
                    f = cs * e[j] + sn * V[j + 1]
                    V[j + 1] = -sn * e[j] + cs * V[j + 1]
                    g = sn * e[j + 1]
                    e[j + 1] = cs * e[j + 1]
                    if (j < m - 1):
                        for i in range(m):
                            t = cs * S[i][j] + sn * S[i][j + 1]
                            S[i][j + 1] = -sn * S[i][j] + cs * S[i][j + 1]
                            S[i][j] = t
                e[p - 2] = f
                itr = itr + 1
            elif caso == 4:
                if V[k] <= 0.0:
                    if V[k] < 0.0:
                        V[k] = -V[k]
                    else:
                        V[k] = 0.0
                    for i in range(pp + 1):
                        D[i][k] = -D[i][k]
                # sort values
                while (k < pp) and (V[k] < V[k + 1]):
                    (V[k], V[k + 1]) = (V[k + 1], V[k])
                    if k < n - 1:
                        for i in range(n):
                            (D[i][k], D[i][k + 1]) = (D[i][k + 1], D[i][k])
                    if k < m - 1:
                        for i in range(m):
                            (S[i][k], S[i][k + 1]) = (S[i][k + 1], S[i][k])
                    k += 1
                itr = 0
                p -= 1

        VM = Matriz.zeros(len(V))
        for l in range(len(V)):
            VM[l][l] = V[l]

        return Matriz(S), VM, Matriz(D)

    def limite_singular(self):
        C = Matriz.zeros(self.m)
        for i in range(self.m):
            C[i][i] = 1/ self[i][i]

        return C

    def solveSVD(self, b):
        (u, s, v) = self.SVD()
        l = s.limite_singular()
        x1 = v * l
        x2 = u.T() * Matriz(b).T()
        x = x1 * x2
        return x.T()[0]

    def faca_zero(self, inplace=True):
        """Converte todos os valores menores que precisao para 0"""

        # Inplace or not?
        if inplace == True:
            M = self
        else:
            M = self.copia()

        m, n = M.m, M.n

        for i in range(m):
            for j in range(n):
                if comozero(M[i][j]):
                    M[i][j] = 0

        return M

    def fatorescala(self):
        """Retorna o vetor com o maximo de cada linha da matriz
        
        s_i = max_j | A_{ij} |, i = 1, 2, 3, \ldots, n

        ref.: [Kiusalaas, 2013] Pag. 71
        """
        m, n = self.m, self.n
        M = self

        s=[]
        for i in range(m):
            s.append(max(abs(M[i])))

        return Vetor(s)

    def tamanhorelativo(self):
        """Retorna a matriz de tamanho relativo

        r_ij = \frac{ | A_ij | }{s_i}

        ref.: [Kiusalaas, 2013] Pag. 71
        """

        M = self.copia()
        m, n = M.m, M.n

        s = M.fatorescala()

        for i in range(m):
            M[i] = abs(M[i]) * (1/s[i])

        return M
        
    def trocalinha(self, i, j, inplace=True):
        """Retorna copia da matriz com linha trocada"""

        if inplace == True:
            M = self
        else:
            M = self.copia()

        m = M.m

        assert 0 <= i < m, "Indice linha deve ser valido"

        M[i], M[j] = M[j], M[i]

        return M

    def trocacoluna(self, i, j, inplace=True):
        """Retorna copia da matriz com coluna trocada"""

        if inplace == True:
            M = self
        else:
            M = self.copia()

        m, n = M.m, M.n

        assert 0 <= i < n, "Indice coluna deve ser valido"

        for l in range(m):
            M[l][i], M[l][j] = M[l][j], M[l][i]

        return M

    def e_simetrica(self):
        """Verifica se a matriz é simetrica"""

        assert(self.e_quadrada), "Matriz nao e quadrada"

        M = self

        for i in range(M.n):
            for j in range(i,M.m):
               if ( M[i][j] != M[j][i] ):
                   return False

        return True

    def normaliza(self):
        """ Retorna matriz normalizada """

        M = self

        for i in range(M.n):
            M[i] = M[i] * (1/M[i].norma())

        return M


"""Classes para tratamento de erros"""


class MatrizError(Exception):
    """Standard Error class for Matriz object"""
    pass


if __name__ == "__main__":
    """Just some initializations to use on iPython tests"""

    # Matriz A para brincar
    A = Matriz([[1, 4, 1], [1, 6, -1], [2, -1, 2]])
    b = Vetor([7,13,5])
    Ab = A.aumenta(b)
    Abg = Ab.gauss()

    # Matriz M para brincar
    M = Matriz([[2, 1, 3], [2, 6, 8], [6, 8, 18]])
    n = Vetor([1, 3, 5])
    Mn = M.aumenta(n)

    # Algumas extracoes de A
    vl = A.linha(0)
    vc = A.coluna(1)
    S = A.submatriz()
     
    # Vandermode 
    vv = Vetor([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    VM = vandermode(vv)
    vb = Vetor([0, 1, 0, 1, 0, 1])
    VMvb = VM.aumenta(vb)


    AW =  Matriz([ 
        [      12,      -51,        4],
        [       6,      167,      -68],
        [      -4,       24,      -41] ])

