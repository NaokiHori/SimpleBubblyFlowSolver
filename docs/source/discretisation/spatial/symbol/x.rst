####################
Wall-normal operator
####################

*******
Average
*******

Average at wall-normal cell faces:

.. math::

   \vat{
      \ave{
         q
      }{
         \gx
      }
   }{
      i + \frac{1}{2}
   }
   =
   \left\{
      \begin{alignedat}{2}
         & \text{Negative wall:} & \vat{q}{\frac{1}{2}}, \\
         & \text{Positive wall:} & \vat{q}{\nx + \frac{1}{2}}, \\
         & \text{Otherwise:} & \frac{1}{2} \vat{q}{i} + \frac{1}{2} \vat{q}{i + 1}
      \end{alignedat}
   \right.

Average at wall-normal cell centers:

.. math::

   \vat{
      \ave{
         q
      }{
         \gx
      }
   }{
      i
   }
   =
   \frac{1}{2} \vat{q}{i - \frac{1}{2}} + \frac{1}{2} \vat{q}{i + \frac{1}{2}}

***************
Differentiation
***************

Differentiation at wall-normal cell faces:

.. math::

   \vat{
      \dif{
         q
      }{
         \gx
      }
   }{
      i + \frac{1}{2}
   }
   =
   \left\{
      \begin{alignedat}{2}
         & \text{Negative wall:} & - \vat{q}{\frac{1}{2}} + \vat{q}{1}, \\
         & \text{Positive wall:} & - \vat{q}{\nx}+ \vat{q}{\nx + \frac{1}{2}} , \\
         & \text{Otherwise:} & - \vat{q}{i} + \vat{q}{i + 1}
      \end{alignedat}
   \right.

Differentiation at wall-normal cell centers:

.. math::

   \vat{
      \dif{
         q
      }{
         \gx
      }
   }{
      i
   }
   =
   - \vat{q}{i - \frac{1}{2}}
   + \vat{q}{i + \frac{1}{2}}

*********
Summation
*********

Summation of a quantity defined at wall-normal cell faces:

.. math::

   \sumxf q
   \equiv
   \vat{q}{\frac{1}{2}}
   +
   \vat{q}{\frac{3}{2}}
   +
   \cdots
   +
   \vat{q}{\nx - \frac{1}{2}}
   +
   \vat{q}{\nx + \frac{1}{2}}

Summation of a quantity defined at wall-normal cell centers:

.. math::

   \sumxc q
   \equiv
   \vat{q}{1}
   +
   \vat{q}{2}
   +
   \cdots
   +
   \vat{q}{\nx - 1}
   +
   \vat{q}{\nx}

