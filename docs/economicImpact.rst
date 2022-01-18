Economic impact model
======================

Data
-----------------------------
* IO Tables
* Table 2 of `Garrett-Peltier, H. (2017) <https://www.sciencedirect.com/science/article/abs/pii/S026499931630709X>`_
* Investment tables (Investment shock per sector and per country, per scenario)

Model
-----------------------------
New Industry: The Final-Demand Approach (Miller y Blair, 2009, cap 13)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Consideremos una economía con dos sectores:

.. math::

   A=\begin{bmatrix}
      a_{11} & a_{12}\\
      a_{21} & a_{22}
   \end{bmatrix}

Pensemos que se incorpora una nueva industria (sector 3). Asumamos que es posible estimar los insumos de los sectores 1 y 2 del valor de la producción del nuevo sector 3; esto es, :math:`a_{13}` y :math:`a_{23}`.

Para cuantificar el impacto de la entrada del sector 3 en la economía, debemos tener alguna medida de la **magnitud** de la nueva actividad económica asociada con el sector 3.

En términos de IO, esto significa que debemos especificar:

* El nivel de producción (producto bruto) del sector 3, :math:`x_3`, o
* La demanda final, :math:`f_3`

Para el ejemplo, se define el nivel de producción de sector 3, que se denota como :math:`\bar{x}_3`.


La nueva demanda de los sectores 1 y 2 que surge por la producción del nuevo sector 3 es :math:`a_{13}\bar{x}_3`  y :math:`a_{23}\bar{x}_3` , respectivamente.

Es decir, podemos ver esas nuevas demandas como un cambio **exógeno** impuesto a los dos sectores originales;

.. math::
  \Delta \mathbf{f}= \begin{bmatrix}
  a_{13}\bar{x}_3 \\
  a_{23}\bar{x}_3
  \end{bmatrix}


de manera que los impactos, en términos del producto de esos dos sectores, estarían dados por :math:`\Delta\mathbf{x} = \mathbf{L}\Delta \mathbf{f}`:

.. math::
  \Delta\mathbf{x}=\begin{bmatrix}
  l_{11} & l_{12}\\
  l_{21} & l_{22}
  \end{bmatrix}
  \begin{bmatrix}
  a_{13}\bar{x}_3 \\
  a_{23}\bar{x}_3
  \end{bmatrix}
  = \begin{bmatrix}
  l_{11}a_{13}\bar{x}_3 + l_{12}a_{23}\bar{x}_3 \\
  l_{21}a_{13}\bar{x}_3 + l_{22}a_{23}\bar{x}_3 \\
  \end{bmatrix}



Dado que también hay una demanda usual, independiente de la demanda del nuevo sector 3, :math:`\bar{f}_1` y :math:`\bar{f}_2`, para esos dos sectores,
el producto bruto total en los dos sectores sería

.. math::
  \begin{bmatrix}
  x_1 \\
  x_2
  \end{bmatrix} =\begin{bmatrix}
  l_{11} & l_{12}\\
  l_{21} & l_{22}
  \end{bmatrix}
  \begin{bmatrix}
  \bar{f}_1 + a_{13}\bar{x}_3 \\
  \bar{f}_2 + a_{23}\bar{x}_3
  \end{bmatrix}
  = \begin{bmatrix}
  l_{11}(\bar{f}_1 + a_{13}\bar{x}_3) + l_{12}(\bar{f}_2 + a_{23}\bar{x}_3) \\
  l_{21}(\bar{f}_1+a_{13}\bar{x}_3) + l_{22}(\bar{f}_2 + a_{23}\bar{x}_3) \\
  \end{bmatrix}


cuando :math:`\bar{f}_1 = 0` y :math:`\bar{f}_2=0`, aislamos el impacto de incoporar el nuevo sector.

Ejemplo 1
"""""""""""

Sea

.. math::
  \mathbf{A}=\begin{bmatrix}
      0.15 & 0.25\\
      0.20& 0.05
  \end{bmatrix},

entonces :math:`(\mathbf{I} - \mathbf{A})^{-1}` es igual a:

.. math::
  \mathbf{A}=\begin{bmatrix}
      1.25412541 & 0.330033\\
      0.2640264  & 1.12211221
  \end{bmatrix}

::

  import numpy as np

  A= np.array([[0.15,0.25],[0.20,0.05]])
  L = np.linalg.inv(np.identity(2)-A)
  L
  >> array([[1.25412541, 0.330033  ],[0.2640264 , 1.12211221]])

Asumamos que las estimaciones de las estimaciones de insumos directos del sector 3 son:

* :math:`a_{13}=0.30`
* :math:`a_{23} = 0.18`

y que se espera que el sector 3 produzca un nivel de 100,000 por año.

De manera que :math:`\bar{x}_3 = 100000`

.. math::
  \Delta \mathbf{f}= \begin{bmatrix}
  0.30 \times 100000  \\
  0.18 \times 100000
  \end{bmatrix}
  =
   \begin{bmatrix}
  30000  \\
  18000
  \end{bmatrix}


El impacto de incoporar el nuevo sector es igual a

.. math::
   \Delta \mathbf{x} = \begin{bmatrix}
  43564  \\
  28118
  \end{bmatrix}

::

  x_bar_3 = 100000
  delta_f = np.array([x_bar_3 * 0.30 ,x_bar_3 * 0.18])
  L@delta_f
  >> array([43564.35643564, 28118.81188119])

El sector 1, al satisfacer la nueva demanda de su producto por un valor de 30,000,
finalmente tendrá que aumentar su producción en 43,560. De manera similar, las nuevas
demandas del sector 2 del sector 3 son de 18,000, pero al final el sector 2 necesitará
producir un total de 28,116 más de producción. Estas cifras representan una forma de
medir el impacto en una economía que surge del movimiento de nueva actividad industrial.

Impactos en el empleo
"""""""""""""""""""""""

Sea :math:`\mathscr{E}` el empleo total y :math:`E = [e_1,e_2,\dots,e_n]` un vector fila
de los coeficientes de trabajo o razones empleo/producto(bruto) de cada sector, la expresión de empleo total sería:

.. math::
  \begin{equation}
  \mathscr{E} = EX
  \end{equation}

Suponiendo una economía con dos sectores, el **impacto** (cambio directo más el indirecto) en el empleo por el incremento exógeno en la demanda final del sector 2 sería igual a

.. math::
  \Delta\mathscr{E}_{d} =
  \begin{bmatrix}
      e_1 & e_2\\
  \end{bmatrix}
  \begin{bmatrix}
      l_{11} & l_{12}\\
      l_{12} & l_{22}
  \end{bmatrix}
  \begin{bmatrix}
      \Delta f_{1} \\
      \Delta f_{2}
  \end{bmatrix}
  =
  E (\mathbf{I} - \mathbf{A})^{-1} \Delta \mathbf{f}
  = E \Delta X

El cambio directo en el empleo debido al incremento en la demanda es :math:`\Delta\mathscr{E}_{d'}`

.. math::
  \Delta\mathscr{E}_{d'} =
  \begin{bmatrix}
      e_1 & e_2\\
  \end{bmatrix}
  \begin{bmatrix}
      \Delta f_{1} \\
      \Delta f_{2}
  \end{bmatrix}
  =
  E\Delta \mathbf{f}



Ejemplo 2
"""""""""""""

Continuando con el ejemplo anterior, calculamos el impacto en el empleo debido a la incorporación del nuevo sector.

Los  coeficientes de trabajo están dados por:

.. math::
  E = \begin{bmatrix}
      0.25 & 0.15\\
  \end{bmatrix}

El cambio en el empleo igual al trabajo requerido por el cambio en la demanda final; es decir

.. math::
  \Delta\mathscr{E}_{d} =
  E \Delta X
  =
  E (\mathbf{I} - \mathbf{A})^{-1} \Delta \mathbf{f}

.. math::
  \Delta\mathscr{E}_{d} =
  \begin{bmatrix}
  0.25 & 0.15
  \end{bmatrix}
  \begin{bmatrix}
  43564  \\
  28118
  \end{bmatrix}
  =
  15108

::

  E = np.array([0.25,0.15])
  E.dot(L@delta_f)
  >> 15108.910891089108

Para calcular el cambio porcentual en el empleo, tenemos que calcular :math:`X`. Suponemos que :math:`\bar{f}_1 = 120000` y :math:`\bar{f}_2=90000`.
Recordando la expresión para calcular :math:`X`;

.. math::
  \begin{bmatrix}
  x_1 \\
  x_2
  \end{bmatrix} =\begin{bmatrix}
  l_{11} & l_{12}\\
  l_{21} & l_{22}
  \end{bmatrix}
  \begin{bmatrix}
  \bar{f}_1  \\
  \bar{f}_2
  \end{bmatrix}

::

  x_bar_3 = 100000
  f1_usual = 120000
  f2_usual = 90000
  X = np.array([f1_usual ,f2_usual])
  (E.dot(L@delta_f)/sum(L@X))*100
  >> 4.829113924050633
