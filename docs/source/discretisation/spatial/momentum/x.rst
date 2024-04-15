#####################
Wall-normal direction
#####################

.. math::

   \pder{\rho \ux}{t}
   =
   &
   \dmomadv{\vx}{\vx}
   \dmomadv{\vx}{\vy}
   \dmomadv{\vx}{\vz}

   &
   \dmompre{\vx}

   &
   \dmomdif{\vx}{\vx}
   \dmomdif{\vx}{\vy}
   \dmomdif{\vx}{\vz}

*********
Advection
*********

.. math::

   \dmomadv{\vx}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is advected in x

.. math::

   \dmomadv{\vx}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is advected in y

.. math::

   \dmomadv{\vx}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is advected in z

*****************
Pressure-gradient
*****************

.. math::

   \dmompre{\vx}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: pressure-gradient contribution

*********
Diffusion
*********

.. math::

   \dmomdif{\vx}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is diffused in x

.. math::

   \dmomdif{\vx}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is diffused in y

.. math::

   \dmomdif{\vx}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
   :language: c
   :tag: ux is diffused in z

