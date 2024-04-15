
.. include:: /reference/reference.txt

#######################
Temporal discretisation
#######################

:ref:`The mass conservation <discrete_mass>` is solved first by means of the :ref:`THINC method <numerical_method>`, which is followed by integrating :ref:`the momentum balance <discrete_momentum>`.

To enforce :ref:`the incompressibility constraint <discrete_incompressibility>`, the updated velocity field is corrected:

.. math::

   \newcommand{bef}{n-\frac{1}{2}}
   \newcommand{aft}{n+\frac{1}{2}}
   \frac{
      u_i^{n+1}
      -
      u_i^*
   }{\Delta t}
   =
   -
   \frac{1}{\rho^{n+1}}
   \frac{1}{h_{\gcs^i}}
   \pder{\psi^{\aft}}{\gcs^i}

where :math:`\psi` is a scalar potential obtained by solving the variable-coefficient Poisson equation:

.. math::

   \frac{1}{J}
   \pder{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      \frac{1}{\rho^{n+1}}
      \frac{1}{h_{\gcs^i}}
      \pder{\psi^{\aft}}{\gcs^i}
   \right)
   =
   \frac{1}{\Delta t}
   \frac{1}{J}
   \pder{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      u_i^*
   \right),

or discretely:

.. math::

   \frac{1}{J}
   \dif{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      \frac{1}{\rho^{n+1}}
      \frac{1}{h_{\gcs^i}}
      \dif{\psi^{\aft}}{\gcs^i}
   \right)
   =
   \frac{1}{\Delta t}
   \frac{1}{J}
   \dif{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      u_i^*
   \right).

Since this is a variable-coefficient Poisson equation, I adopt the approach proposed by |DODD2014| to utilise the orthogonal decomposition.
With this approximation and :math:`\rho_{ref} \equiv \min \left( 1, \hat{\rho} \right)`, the Poisson equation is modified as

.. math::

   \frac{1}{J}
   \dif{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      \frac{1}{h_{\gcs^i}}
      \dif{\psi^{\aft}}{\gcs^i}
   \right)
   =
   \frac{1}{J}
   \dif{}{\gcs^i}
   \left\{
      \frac{J}{h_{\gcs^i}}
      \left(
         \frac{\rho_{ref}}{\rho^{n+1}}
         -
         1
      \right)
      \frac{1}{h_{\gcs^i}}
      \dif{\psi^{\bef}}{\gcs^i}
   \right\}
   +
   \frac{\rho_{ref}}{\Delta t}
   \frac{1}{J}
   \dif{}{\gcs^i}
   \left(
      \frac{J}{h_{\gcs^i}}
      u_i^*
   \right).

.. myliteralinclude:: /../../src/fluid/compute_potential.c
   :language: c
   :tag: additional contribution

The velocity corrections are modified as

.. math::

   \newcommand{\old}[1]{
      +
      \left(
         \frac{1}{\rho^{n+1}}
         -
         \frac{1}{\rho_{ref}}
      \right)
      \frac{1}{h_{\gcs^#1}}
      \pder{\psi^{\bef}}{\gcs^#1}
   }
   \newcommand{\new}[1]{
      -
      \frac{1}{\rho_{ref}}
      \frac{1}{h_{\gcs^#1}}
      \pder{\psi^{\aft}}{\gcs^#1}
   }
   \frac{
      u_i^{n+1}
      -
      u_i^*
   }{\Delta t}
   =
   \new{i}
   \old{i}.

There are two contributions, which are from the new scalar potential and the old scalar potential.

The contributions of the new scalar potential:

.. math::

   \new{\vx}

.. myliteralinclude:: /../../src/fluid/correct_velocity/ux.c
   :language: c
   :tag: new scalar potential contribution

.. math::

   \new{\vy}

.. myliteralinclude:: /../../src/fluid/correct_velocity/uy.c
   :language: c
   :tag: new scalar potential contribution

.. math::

   \new{\vz}

.. myliteralinclude:: /../../src/fluid/correct_velocity/uz.c
   :language: c
   :tag: new scalar potential contribution

The contributions of the old scalar potential:

.. math::

   \old{\vx}

.. myliteralinclude:: /../../src/fluid/correct_velocity/ux.c
   :language: c
   :tag: old scalar potential contribution

.. math::

   \old{\vy}

.. myliteralinclude:: /../../src/fluid/correct_velocity/uy.c
   :language: c
   :tag: old scalar potential contribution

.. math::

   \old{\vz}

.. myliteralinclude:: /../../src/fluid/correct_velocity/uz.c
   :language: c
   :tag: old scalar potential contribution

Since I treat the diffusive terms fully-explicitly, the scalar potential and the pressure are simply related by

.. math::

   p^{n+1}
   -
   p^n
   =
   \psi^{\aft}.

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: explicit contribution

This is followed by updating :math:`\psi`:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: update psi

