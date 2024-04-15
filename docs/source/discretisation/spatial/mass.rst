
.. _discrete_mass:

#################
Mass conservation
#################

The mass conservation is defined at each cell center and described as:

.. math::

   \pder{\rho}{t}
   +
   \frac{1}{J}
   \dif{}{\gx}
   \left(
      \jhx
      \rho
      \ux
   \right)
   +
   \frac{1}{J}
   \dif{}{\gy}
   \left(
      \jhy
      \rho
      \uy
   \right)
   +
   \frac{1}{J}
   \dif{}{\gz}
   \left(
      \jhz
      \rho
      \uz
   \right)
   =
   0.

Specifically the rescaled version

.. math::

   \pder{H}{t}
   +
   \frac{1}{J}
   \dif{}{\gx}
   \left(
      \jhx
      \ux
      H
   \right)
   +
   \frac{1}{J}
   \dif{}{\gy}
   \left(
      \jhy
      \uy
      H
   \right)
   +
   \frac{1}{J}
   \dif{}{\gz}
   \left(
      \jhz
      \uz
      H
   \right)
   =
   0

is considered.

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: compute source of volume-of-fluid

