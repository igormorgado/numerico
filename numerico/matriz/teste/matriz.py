from numerico.matriz import *

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

