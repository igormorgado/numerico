import numpy as np

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


