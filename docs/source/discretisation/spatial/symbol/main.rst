######
Symbol
######

Averages, differentiations, and summations are denotes as follows.

.. toctree::
   :maxdepth: 1

   x
   y
   z

**********************************************
Discrete scale factor and Jacobian determinant
**********************************************

In Cartesian coordinate systems, scale factors are simply equal to the grid sizes:

.. math::

   \hx
   &
   \equiv
   \Delta x,

   \hy
   &
   \equiv
   \Delta y,

   \hz
   &
   \equiv
   \Delta z.

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: scale factor in x, defined at x cell faces

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: scale factor in x, defined at x cell centers

Note that the wall-normal grid sizes are halved on the walls:

.. math::

   \vat{\Delta x}{\frac{1}{2}}
   =
   \vat{\Delta x}{\nx + \frac{1}{2}}
   =
   \frac{1}{2} \Delta x_{\text{rest}}.

The Jacobian determinants, which are also defined at cell centers and faces, are simply the product of the local scale factors:

.. math::

   J
   \equiv
   \hx
   \hy
   \hz.

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: Jacobian determinant, defined at x cell faces

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: Jacobian determinant, defined at x cell centers

