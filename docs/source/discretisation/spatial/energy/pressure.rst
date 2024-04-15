The pressure-gradient terms

.. math::

   -
   \frac{1}{h_{\gcs^i}}
   \dif{p}{\gcs^i},

which contributes to the energy balance as follows:

.. math::

   -
   \sumzc
   \sumyc
   \sumxf
   J
   \ux
   \frac{1}{\hx}
   \dif{p}{\gx}
   =
   \sumzc
   \sumyc
   \sumxc
   J
   p
   \frac{1}{J}
   \dif{
      \left(
         \jhx
         \ux
      \right)
   }{\gx}

.. math::

   -
   \sumzc
   \sumyf
   \sumxc
   J
   \uy
   \frac{1}{\hy}
   \dif{p}{\gy}
   =
   \sumzc
   \sumyc
   \sumxc
   J
   p
   \frac{1}{J}
   \dif{
      \left(
         \jhy
         \uy
      \right)
   }{\gy}

.. math::

   -
   \sumzf
   \sumyc
   \sumxc
   J
   \uz
   \frac{1}{\hz}
   \dif{p}{\gz}
   =
   \sumzc
   \sumyc
   \sumxc
   J
   p
   \frac{1}{J}
   \dif{
      \left(
         \jhz
         \uz
      \right)
   }{\gz}

The sum is

.. math::

   \sumzc
   \sumyc
   \sumxc
   J
   p
   \left\{
      \frac{1}{J}
      \dif{
         \left(
            \jhx
            \ux
         \right)
      }{\gx}
      +
      \frac{1}{J}
      \dif{
         \left(
            \jhy
            \uy
         \right)
      }{\gy}
      +
      \frac{1}{J}
      \dif{
         \left(
            \jhz
            \uz
         \right)
      }{\gz}
   \right\},

which is zero because the component inside the wavy parentheses is :ref:`the incompressibility constraint <discrete_incompressibility>`.

