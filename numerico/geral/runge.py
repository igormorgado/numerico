import numpy as np

def plota_runge(x, cor=None, legenda=None):
    """Retorna o grafico do Fenomeno de Runge nos pontos"""
    
    
    plotargs = { 'alpha': .5,
                 'marker': 'o',
                 'linestyle': ''
               }
    
    if cor is not None:
        plotargs['color'] = cor

    if legenda is not None:
        plotargs['label'] = legenda
    
    return plt.plot(x, runge(x), **plotargs)


def runge(x):
    """Runge function"""
    return 1 / (1 + 25 * x ** 2)
