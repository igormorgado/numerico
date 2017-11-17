import numpy as np
import geral.runge

def spline_linear(x, f, x_k, x_ki):
    """Returns the linear interpolation of f(x) beteween x_k and x_ki."""
    A = (x_ki - x) / (x_ki - x_k)
    B = (x - x_k)  / (x_ki - x_k)
    
    return A*f(x_k) + B*f(x_ki)

def plota_spline_linear(x, f, pontos=10, legenda=None):
    
    z = np.empty(0)
    l = np.empty(0)
    
    for i in range(len(x)-1):
        z_i = np.linspace(x[i],x[i+1],pontos)[1:-1]
        l_i = spline_linear(z_i, runge, x[i], x[i+1])
        z = np.concatenate((z, z_i))
        l = np.concatenate((l, l_i))

    plotargs = {}
    if legenda is not None:
        plotargs['label'] = legenda
        
    return plt.plot(z, l, marker='.', linestyle='', alpha=.1, **plotargs)    

