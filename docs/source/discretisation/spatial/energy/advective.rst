I consider the effects of

.. math::

   \pder{\rho u_i}{t}
   +
   \frac{1}{J}
   \dif{
      \left(
         \ave{
            \frac{J}{h_{\gcs^j}}
            \rho u_j
         }{\gcs^i}
         \ave{
            u_i
         }{\gcs^j}
      \right)
   }{\gcs^j}

on the global kinetic energy balance; namely:

.. math::

   \sum_{i\,\text{face}}
   J
   u_i
   \left\{
      \pder{\rho u_i}{t}
      +
      \frac{1}{J}
      \dif{
         \left(
            \ave{
               \frac{J}{h_{\gcs^j}}
               \rho u_j
            }{\gcs^i}
            \ave{
               u_i
            }{\gcs^j}
         \right)
      }{\gcs^j}
   \right\}.

Contribution of the first term leads to

.. math::

   \sum_{i\,\text{face}}
   J
   u_i
   \pder{\rho u_i}{t}
   =
   \sum_{i\,\text{face}}
   J
   \pder{\rho u_i u_i}{t}
   -
   \sum_{i\,\text{face}}
   J
   \rho u_i
   \pder{u_i}{t},

while the second term gives

.. math::

   \sum_{i\,\text{face}}
   J
   u_i
   \frac{1}{J}
   \dif{
      \left(
         \ave{
            \frac{J}{h_{\gcs^j}}
            \rho u_j
         }{\gcs^i}
         \ave{
            u_i
         }{\gcs^j}
      \right)
   }{\gcs^j}
   &
   =
   -
   \sum_{i\,\text{center}}
   \ave{
      u_i
   }{\gcs^j}
   \ave{
      \frac{J}{h_{\gcs^j}}
      \rho u_j
   }{\gcs^i}
   \dif{
      u_i
   }{\gcs^j} \\
   &
   =
   -
   \sum_{i\,\text{face}}
   J
   u_i
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\gcs^j}}
         \rho u_j
      }{\gcs^i}
      \dif{
         u_i
      }{\gcs^j}
   }{\gcs^j}.

As a result, I obtain

.. math::

   \sum_{i\,\text{face}}
   J
   \pder{\rho u_i u_i}{t}
   -
   \sum_{i\,\text{face}}
   J
   u_i
   \left(
      \rho
      \pder{u_i}{t}
      +
      \frac{1}{J}
      \ave{
         \ave{
            \frac{J}{h_{\gcs^j}}
            \rho u_j
         }{\gcs^i}
         \dif{
            u_i
         }{\gcs^j}
      }{\gcs^j}
   \right).

What I have inside the parentheses of the second term is the so-called gradient form of the advective terms.

