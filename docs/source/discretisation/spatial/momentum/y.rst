#####################
Stream-wise direction
#####################

.. math::

   \pder{\rho \uy}{t}
   =
   &
   \dmomadv{\vy}{\vx}
   \dmomadv{\vy}{\vy}
   \dmomadv{\vy}{\vz}

   &
   \dmompre{\vy}

   &
   \dmomdif{\vy}{\vx}
   \dmomdif{\vy}{\vy}
   \dmomdif{\vy}{\vz}

*********
Advection
*********

.. math::

   \dmomadv{\vy}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is advected in x

.. math::

   \dmomadv{\vy}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is advected in y

.. math::

   \dmomadv{\vy}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is advected in z

*****************
Pressure-gradient
*****************

.. math::

   \dmompre{\vy}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: pressure-gradient contribution

*********
Diffusion
*********

.. math::

   \dmomdif{\vy}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is diffused in x

.. math::

   \dmomdif{\vy}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is diffused in y

.. math::

   \dmomdif{\vy}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
   :language: c
   :tag: uy is diffused in z

