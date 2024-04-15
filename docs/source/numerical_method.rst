
.. _numerical_method:

.. include:: /reference/reference.txt

################
Numerical method
################

*******************
Interface capturing
*******************

To begin with, I consider an indicator function :math:`H`, which works as a marker taking :math:`0` or :math:`1` when the infinitesimal control volume is occupied with the primary or the secondary liquids, respectively.
By this definition, as an equation governing the evolution of :math:`H`, I have

.. math::

   \pder{H}{t}
   +
   \frac{u_j}{h_{\gcs^j}}
   \pder{H}{\gcs^j}
   =
   0.

Since I assume the liquids are incompressible, this advection equation is equal to

.. math::

   \pder{H}{t}
   +
   \frac{1}{J}
   \pder{}{\gcs^j}
   \left(
      \frac{J}{h_{\gcs^j}}
      u_j
      H
   \right)
   =
   0,

whose volume-integrated form leads to

.. math::

   J
   \pder{\phi}{t}
   +
   \int_{\partial V}
   \frac{J}{h_{\gcs^j}}
   u_j
   H
   n_j
   dS
   =
   0,

where :math:`\phi` is the rate of the primary phase volume inside the control volume, or known as `volume-of-fluid`.
Note that I assume the coordinate systems do not change in time.
To integrate this equation, I employ the THINC/QQ method (|XIE2017|), whose numerical treatment is extensively discussed in `the other project <https://naokihori.github.io/SimpleVOFSolver/numerical_method/main.html>`_.

***********
Consistency
***********

Note that the following idea is widely known in the phase-field community (e.g. |MIRJALILI2021|), which is slightly customised and applied here.

As discussed above, how the volume fraction is evolved is totally determined by the volume-of-fluid method.
Since :math:`H` and :math:`\rho` are related by

.. math::

   \rho = 1 + \left( \hat{\rho} - 1 \right) H,

the evolution of the density (i.e. mass per unit volume) is also governed by the VOF method.
Since the momentum is tightly linked to the flux of the mass, it is natural to apply this information (mass flux) to update the momentum field.
To this end, I focus on the evolution of :math:`\phi` again:

.. math::

   J
   \pder{\phi}{t}
   +
   \pder{}{\gcs^j}
   \left(
      \frac{J}{h_{\gcs^j}}
      H
      u_j
   \right)
   =
   0.

As the flux of the indicator function :math:`H u_j` is already known, the density flux :math:`\rho u_j` is obtained by

.. math::

   \rho u_j
   =
   u_j
   +
   \left(
      \hat{\rho}
      -
      1
   \right)
   H
   u_j.

.. myliteralinclude:: /../../src/interface/mass_flux.c
   :language: c
   :tag: convert x vof flux to x mass flux

.. myliteralinclude:: /../../src/interface/mass_flux.c
   :language: c
   :tag: convert y vof flux to y mass flux

.. myliteralinclude:: /../../src/interface/mass_flux.c
   :language: c
   :tag: convert z vof flux to z mass flux

*********************
Surface tension force
*********************

For simplicity, the continuum surface force model (|BRACKBILL1992|) is adopted to model the interfacial tension:

.. math::

   f_i
   \approx
   \frac{2 \rho}{1 + \hat{\rho}}
   \sigma
   \kappa
   \frac{1}{h_{\gcs^i}}
   \pder{\phi}{\gcs^i}.

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in x direction

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in y direction

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in z direction

The pre-factor is computed here:

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute density factor

