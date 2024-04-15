###################
Span-wise direction
###################

.. math::

   \pder{\rho \uz}{t}
   =
   &
   \dmomadv{\vz}{\vx}
   \dmomadv{\vz}{\vy}
   \dmomadv{\vz}{\vz}

   &
   \dmompre{\vz}

   &
   \dmomdif{\vz}{\vx}
   \dmomdif{\vz}{\vy}
   \dmomdif{\vz}{\vz}

*********
Advection
*********

.. math::

   \dmomadv{\vz}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is advected in x

.. math::

   \dmomadv{\vz}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is advected in y

.. math::

   \dmomadv{\vz}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is advected in z

*****************
Pressure-gradient
*****************

.. math::

   \dmompre{\vz}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: pressure-gradient contribution

*********
Diffusion
*********

.. math::

   \dmomdif{\vz}{\vx}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is diffused in x

.. math::

   \dmomdif{\vz}{\vy}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is diffused in y

.. math::

   \dmomdif{\vz}{\vz}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
   :language: c
   :tag: uz is diffused in z

