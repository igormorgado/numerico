import numpy as np
import diferencas.divdiff

def fwd_diff(y):
    """Returns a forward difference of a set of points"""

    dd = np.zeros([len(y), len(y)]) 
    dd[0] = y
    dd[1,:-1] = (dd[0,1:] - dd[0,:-1])

    for i in range(2,len(y)):
        dd[i,:-i] = dd[i-1,1:-(i-1)] - dd[i-1,:-i]

    return dd


def fwd_div_diff(x,y):
    """Returns the forward divided difference of set of points"""

    dd = np.zeros([len(y), len(y)]) 
    dd[0] = y
    dd[1,:-1] = (dd[0,1:] - dd[0,:-1]) / (x[1:]-x[:-1])

    for i in range(2,len(y)):
        dd[i,:-i] = (dd[i-1,1:-(i-1)] - dd[i-1,:-i]) / (x[i:]-x[:-i])

    return dd


def interp_coeff(x,y):
    """Returns the highest order interpolation coefficients 
    to a set of data points"""

    return fwd_div_diff(x,y)[:,0]

def interpolate(t, x, y):
    """Returns the points 't' interpolated by 'x,y' points.
    The interpolation order will be based on size of 'x,y'"""

    b = interp_coeff(x, y)

    f = np.full_like(t, b[0], dtype='float64')

    for i in range(1,len(b)):
        mul = 1
        for m in range(i):
            mul *= t - x[m]

        f += b[i] * mul

    return f


if __name__ == "__main__":
    x = np.array([0, 10, 15, 20, 22.5, 30])
    y = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.36])

    vx = x[1:4]
    vy = y[1:4]
    t=np.linspace(vx[0],vx[-1],21)
    ft = interpolate(t,vx,vy)
    coeff2 = interp_coeff(vx,vy)

