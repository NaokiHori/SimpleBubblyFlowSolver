
.. _equation:

.. include:: /reference/reference.txt

########
Equation
########

******************
Boundary condition
******************

Note that the following set of boundary conditions is assumed throughout the project.

* The domain is wall-bounded in the :math:`x` direction, while periodic in the other directions.
* The walls are no-slip and impermeable.
* The walls may move in the :math:`y` direction with constant speeds over time, while there is no :math:`z` motion.

******
Metric
******

For simplicity and generality, the governing equations are written in a general rectilinear coordinate system :math:`\gcs^i` with normalised Cartesian components :math:`u_i`.
:math:`h_{\gcs^i}` denote scale factors and its product is the Jacobian determinant :math:`J` due to the orthogonality.

See `e.g. this <https://naokihori.github.io/OrthogonalNS/index.html>`_ for the derivations.

**********************
Non-dimensionalisation
**********************

To normalise the equations, I use the density :math:`\rho` and the dynamic viscosity :math:`\mu` of the primary liquid as reference values.
Namely, the equations described below all assume that the density and the dynamic viscosity of the primary liquids are unity.
In particular, through the phase indicator function :math:`H` (defined later), the local and the instantaneous density and dynamic viscosity are given by

.. math::

   \rho \left( x_i, t \right) & = 1 + \left( \hat{\rho} - 1 \right) H \left( x_i, t \right),

   \mu \left( x_i, t \right) & = 1 + \left( \hat{\mu} - 1 \right) H \left( x_i, t \right).

Additionally, by using the reference length / velocity scales, we have three non-dimensional numbers: :math:`Re, We, Fr`.

.. note::

    To achieve the shear-stress continuity, the dynamic viscosity should be given as

    .. math::

        \frac{1}{\mu \left( x_i, t \right)} = 1 + \left( \frac{1}{\hat{\mu}} - 1 \right) H \left( x_i, t \right),

    as pointed out by e.g., |FERZIGER2003|.
    For the time being, however, I use the arithmetic average to be consistent with the density treatment, which could be controversial.

    I would like to acknowledge Dr. Kevin Zhong for fruitful discussion regarding the viscosity formulation.

************************
Velocity-gradient tensor
************************

.. include:: velocity_gradient_tensor.rst

*******************
Shear-stress tensor
*******************

The shear-stress tensor for Newtonian liquids is defined as

.. math::

   \tau_{ij}
   \equiv
   2 \mu s_{ij},

where :math:`s_{ij}` is the strain-rate tensor:

.. math::

   s_{ij}
   \equiv
   \frac{1}{2}
   l_{ij}
   +
   \frac{1}{2}
   l_{ji}.

Thus

.. math::

   \tau_{ij}
   =
   \mu
   l_{ij}
   +
   \mu
   l_{ji}.

****************************
Incompressibility constraint
****************************

.. include:: incompressibility.rst

*****************
Mass conservation
*****************

.. include:: mass.rst

****************
Momentum balance
****************

.. include:: mom.rst

******************
Quadratic quantity
******************

I consider the quadratic quantities

.. math::

   k_i
   \equiv
   \frac{1}{2}
   \rho
   u_i u_i \,\, \text{(No summation)},

which satisfy the following relations.

.. include:: quad.rst

By volume-integrating these three relations inside the whole domain and summing them up, I obtain the relation of the global kinetic energy:

.. math::

   \pder{}{t}
   \int_V
   \left(
      \kx
      +
      \ky
      +
      \kz
   \right)
   J
   d\gx
   d\gy
   d\gz
   =
   \left( \text{transport} \right)
   +
   \left( \text{dissipation} \right).

Note that the body force contributions :math:`f_i u_i` are omitted.

Here the *transport* is the net kinetic energy going through the walls which attributes to the first diffusive term in the stream-wise momentum equation:

.. math::

   -
   \int_{\gz}
   \int_{\gy}
   \vat{
      \left(
         \jhx
         \uy
         \tau_{\vx \vy}
      \right)
   }{x = 0}
   d\gy
   d\gz
   +
   \int_{\gz}
   \int_{\gy}
   \vat{
      \left(
         \jhx
         \uy
         \tau_{\vx \vy}
      \right)
   }{x = 1}
   d\gy
   d\gz
   =
   -
   \int_{S, x = 0}
   \uy
   \tau_{\vx \vy}
   dS
   +
   \int_{S, x = 1}
   \uy
   \tau_{\vx \vy}
   dS,

while the *dissipation* is handled by the other terms

.. math::

   -
   \int_V
   \left(
      \begin{aligned}
         &
         +
         \frac{1}{\hx}
         \pder{\ux}{\gx}
         \tau_{\vx \vx}
         +
         \frac{1}{\hy}
         \pder{\ux}{\gy}
         \tau_{\vy \vx}
         +
         \frac{1}{\hz}
         \pder{\ux}{\gz}
         \tau_{\vz \vx} \\
         &
         +
         \frac{1}{\hx}
         \pder{\uy}{\gx}
         \tau_{\vx \vy}
         +
         \frac{1}{\hy}
         \pder{\uy}{\gy}
         \tau_{\vy \vy}
         +
         \frac{1}{\hz}
         \pder{\uy}{\gz}
         \tau_{\vz \vy} \\
         &
         +
         \frac{1}{\hx}
         \pder{\uz}{\gx}
         \tau_{\vx \vz}
         +
         \frac{1}{\hy}
         \pder{\uz}{\gy}
         \tau_{\vy \vz}
         +
         \frac{1}{\hz}
         \pder{\uz}{\gz}
         \tau_{\vz \vz}
      \end{aligned}
   \right)
   dV
   =
   -
   \int_V
   \left(
      \begin{aligned}
         &
         +
         l_{\vx \vx}
         \tau_{\vx \vx}
         +
         l_{\vy \vx}
         \tau_{\vy \vx}
         +
         l_{\vz \vx}
         \tau_{\vz \vx} \\
         &
         +
         l_{\vx \vy}
         \tau_{\vx \vy}
         +
         l_{\vy \vy}
         \tau_{\vy \vy}
         +
         l_{\vz \vy}
         \tau_{\vz \vy} \\
         &
         +
         l_{\vx \vz}
         \tau_{\vx \vz}
         +
         l_{\vy \vz}
         \tau_{\vy \vz}
         +
         l_{\vz \vz}
         \tau_{\vz \vz}
      \end{aligned}
   \right)
   dV.

Note that the advective contributions on the global energy balance vanish due to the prescribed boundary conditions.

