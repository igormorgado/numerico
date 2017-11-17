#!/usr/bin/env python3

import numpy as np
from numpy import sum, abs, diag, max, zeros_like
from numpy.linalg import matrix_rank, solve, inv
from erros import eps



def dominancia(A):
    """Retorna o fator de dominancia de cada linha da matriz A"""
    return (sum(abs(A), axis=1) - abs(diag(A))) / abs(diag(A))


if __name__ == "__main__":
    A = np.array([[3, -1, 0,  1],
                  [1,  2, 0,  0],
                  [2,  3, 7, -1],
                  [1, -1, 1,  6]])
    
    D = diag(diag(A))
    B = solve(D, A)
    B2 = inv(D) @ A
    I = eye(matrix_rank(A))
    E = I - B
    
    G = zeros_like(E)
    for i in range(100):
        G += matrix_power(E, i)
    
    R = I - B @ G
    print(all(R < eps))
    
