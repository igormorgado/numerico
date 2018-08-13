#!/usr/bin/env python3

import numpy as np
from numpy import sum, abs, diag 



def dominancia(A):
    """Retorna o fator de dominancia de cada linha da matriz A"""
    return (sum(abs(A), axis=1) - abs(diag(A))) / abs(diag(A))

