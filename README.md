NUMERICO
========

Biblioteca python que implementa metodos numericos comumente vistos em aulas de graduacao e posgraduacao.

Esta biblioteca nao tem como objetivo ser eficiente, rapida ou abrangente, mas sim simples de ser lida e 
utilizada em casos mais gerais.

Tambem pode ser utilizada para estudar os metodos em sala de aula.

Ainda esta em estagio completamente BETA

Para bibliotecas eficientes veja:  numpy, scipy


INSTALAR
--------

```
git clone https://github.com/igormorgado/numerico.git numerico

pip install numerico
```


EXEMPLO
-------

```
from numerico.alglin.vetor import Vetor

u=Vetor([1,2,3])

v=Vetor([4,5,6])

u+v
OUT: numerico.alglin.vetor.Vetor([    5,     7,     9])
```


Jan 2017
