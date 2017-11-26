import numpy as np
from numerico.diferencas.dif_divididas import dif_div_prog
from numerico.interpolacao.newton import *


if __name__ == "__main__":
    x = np.array([0, 10, 15, 20, 22.5, 30])
    y = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.36])

    vx = x[1:4]
    vy = y[1:4]
    t=np.linspace(vx[0],vx[-1],21)
    ft = interpolate(t,vx,vy)
    coeff2 = interp_coeff(vx,vy)

