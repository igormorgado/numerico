from setuptools import setup, find_packages

setup(
    name = 'numerico',
    packages = find_packages(), 
    version = '0.0.1',
    install_requires = [ 'numpy', 'bitstring' ],
    description = 'Metodos numericos para todos',
    url = 'http://github.com/igormorgado/numerico',
    license = 'GNU/GPLv2',
    python_requires = '>=3.4, <4',
    author = 'Igor Morgado',
    author_email = 'morgado.igor@gmail.com',
    keywords = ['latex', 'numpy', 'numerical'],
    classifiers = [
          "Programming Language :: Python",
          "Programming Language :: Python :: 3",
          "Development Status :: 2 - Pre-Alpha",
          "Environment :: Other Environment",
          "Intended Audience :: Science/Research",
          "Intended Audience :: Education",
          "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
          "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    long_description = """\
Ferramentas de calculo numerico de facil acesso
-------------------------------------------------

O objetivo e' criar bibliotecas de calculo numerico com codigo simples 
e acessivel para graduandos lusofonos.
"""

)
