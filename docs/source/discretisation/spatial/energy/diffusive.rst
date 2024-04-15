In general, the diffusive terms are written as

.. math::

   \frac{1}{J}
   \dif{}{\gcs^j}
   \left(
      \frac{J}{h_{\gcs^j}}
      \tau_{i j}
   \right),

which is the diffusion of the :math:`i`-th momentum in the :math:`j`-th direction.

The global contribution

.. math::

   \sum_{i\,\text{face}}
   J
   u_i
   \frac{1}{J}
   \dif{}{\gcs^j}
   \left(
      \frac{J}{h_{\gcs^j}}
      \tau_{i j}
   \right)
   =
   \sum_{i\,\text{face}}
   u_i
   \dif{}{\gcs^j}
   \left(
      \frac{J}{h_{\gcs^j}}
      \tau_{i j}
   \right)

is as follows.

:math:`x` momentum contribution:

.. math::

   &
   \sumzc
   \sumyc
   \sumxf
   \ux
   \dif{}{\gx}
   \left(
      \jhx
      \tau_{\vx \vx}
   \right)
   -
   =
   \sumzc
   \sumyc
   \sumxc
   J
   l_{\vx \vx}
   \tau_{\vx \vx}

   &
   \sumzc
   \sumyc
   \sumxf
   \ux
   \dif{}{\gy}
   \left(
      \jhy
      \tau_{\vx \vy}
   \right)
   =
   -
   \sumzc
   \sumyf
   \sumxf
   J
   l_{\vx \vy}
   \tau_{\vx \vy}

   &
   \sumzc
   \sumyc
   \sumxf
   \ux
   \dif{}{\gz}
   \left(
      \jhz
      \tau_{\vx \vz}
   \right)
   =
   -
   \sumzf
   \sumyc
   \sumxf
   J
   l_{\vx \vz}
   \tau_{\vx \vz}

:math:`y` momentum contribution:

.. math::

   &
   \sumzc
   \sumyf
   \sumxc
   \uy
   \dif{}{\gx}
   \left(
      \jhx
      \tau_{\vy \vx}
   \right)
   =
   -
   \sumzc
   \sumyf
   \vat{
      \left(
         \uy
         \jhx
         \tau_{\vy \vx}
      \right)
   }{\frac{1}{2}}
   +
   \sumzc
   \sumyf
   \vat{
      \left(
         \uy
         \jhx
         \tau_{\vy \vx}
      \right)
   }{\nx + \frac{1}{2}}
   -
   \sumzc
   \sumyf
   \sumxf
   J
   l_{\vy \vx}
   \tau_{\vy \vx}

   &
   \sumzc
   \sumyf
   \sumxc
   \uy
   \dif{}{\gy}
   \left(
      \jhy
      \tau_{\vy \vy}
   \right)
   =
   -
   \sumzc
   \sumyc
   \sumxc
   J
   l_{\vy \vy}
   \tau_{\vy \vy}

   &
   \sumzc
   \sumyf
   \sumxc
   \uy
   \dif{}{\gz}
   \left(
      \jhz
      \tau_{\vy \vz}
   \right)
   =
   -
   \sumzf
   \sumyf
   \sumxc
   J
   l_{\vy \vz}
   \tau_{\vy \vz}

:math:`z` momentum contribution:

.. math::

   &
   \sumzf
   \sumyc
   \sumxc
   \uz
   \dif{}{\gx}
   \left(
      \jhx
      \tau_{\vz \vx}
   \right)
   =
   -
   \sumzf
   \sumyc
   \sumxf
   J
   l_{\vz \vx}
   \tau_{\vz \vx}

   &
   \sumzf
   \sumyc
   \sumxc
   \uz
   \dif{}{\gy}
   \left(
      \jhy
      \tau_{\vz \vy}
   \right)
   =
   -
   \sumzf
   \sumyf
   \sumxc
   J
   l_{\vz \vy}
   \tau_{\vz \vy}

   &
   \sumzf
   \sumyc
   \sumxc
   \uz
   \dif{}{\gz}
   \left(
      \jhz
      \tau_{\vz \vz}
   \right)
   =
   -
   \sumzc
   \sumyc
   \sumxc
   J
   l_{\vz \vz}
   \tau_{\vz \vz}

All terms are dissipative, except the two terms in the :math:`y` contribution which are the energy throughput on the walls.

