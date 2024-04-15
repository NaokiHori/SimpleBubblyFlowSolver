
.. _discrete_incompressibility:

############################
Incompressibility constraint
############################

The incompressibility constraint is defined at each cell center:

.. math::

   \frac{1}{J}
   \dif{}{\gx}
   \left(
      \jhx
      \ux
   \right)
   +
   \frac{1}{J}
   \dif{}{\gy}
   \left(
      \jhy
      \uy
   \right)
   +
   \frac{1}{J}
   \dif{}{\gz}
   \left(
      \jhz
      \uz
   \right)
   =
   0.

.. myliteralinclude:: /../../src/logging/divergence.c
   :language: c
   :tag: compute local divergence

.. myliteralinclude:: /../../src/fluid/compute_potential.c
   :language: c
   :tag: compute local divergence

