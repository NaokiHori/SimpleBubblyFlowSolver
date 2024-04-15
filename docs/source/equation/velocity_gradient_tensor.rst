The velocity-gradient tensor is defined as

.. math::

   \vec{e}_i
   \otimes
   \vec{e}_j
   l_{ij},

where the components are given as

.. math::

   \begin{pmatrix}
      l_{\vx \vx} & l_{\vx \vy} & l_{\vx \vz} \\
      l_{\vy \vx} & l_{\vy \vy} & l_{\vy \vz} \\
      l_{\vz \vx} & l_{\vz \vy} & l_{\vz \vz} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \frac{1}{\hx}
      \pder{\ux}{\gx}
      &
      \frac{1}{\hx}
      \pder{\uy}{\gx}
      &
      \frac{1}{\hx}
      \pder{\uz}{\gx}
      \\
      \frac{1}{\hy}
      \pder{\ux}{\gy}
      &
      \frac{1}{\hy}
      \pder{\uy}{\gy}
      &
      \frac{1}{\hy}
      \pder{\uz}{\gy}
      \\
      \frac{1}{\hz}
      \pder{\ux}{\gz}
      &
      \frac{1}{\hz}
      \pder{\uy}{\gz}
      &
      \frac{1}{\hz}
      \pder{\uz}{\gz}
   \end{pmatrix}.

