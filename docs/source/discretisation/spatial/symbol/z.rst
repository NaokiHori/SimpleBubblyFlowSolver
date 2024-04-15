##################
Span-wise operator
##################

*******
Average
*******

Average at span-wise cell faces:

.. math::

   \vat{
      \ave{
         q
      }{
         \gz
      }
   }{
      k + \frac{1}{2}
   }
   =
   \frac{1}{2} \vat{q}{k} + \frac{1}{2} \vat{q}{k + 1}

Average at span-wise cell centers:

.. math::

   \vat{
      \ave{
         q
      }{
         \gz
      }
   }{
      k
   }
   =
   \frac{1}{2} \vat{q}{k - \frac{1}{2}} + \frac{1}{2} \vat{q}{k + \frac{1}{2}}

***************
Differentiation
***************

Differentiation at span-wise cell faces:

.. math::

   \vat{
      \dif{
         q
      }{
         \gz
      }
   }{
      k + \frac{1}{2}
   }
   =
   -
   \vat{q}{k}
   +
   \vat{q}{k + 1}

Differentiation at span-wise cell centers:

.. math::

   \vat{
      \dif{
         q
      }{
         \gz
      }
   }{
      k
   }
   =
   -
   \vat{q}{k - \frac{1}{2}}
   +
   \vat{q}{k + \frac{1}{2}}

*********
Summation
*********

Summation of a quantity defined at span-wise cell faces:

.. math::

   \sumzf q
   \equiv
   \vat{q}{\frac{1}{2}}
   +
   \vat{q}{\frac{3}{2}}
   +
   \cdots
   +
   \vat{q}{\nz - \frac{3}{2}}
   +
   \vat{q}{\nz - \frac{1}{2}}

Summation of a quantity defined at span-wise cell centers:

.. math::

   \sumzc q
   \equiv
   \vat{q}{1}
   +
   \vat{q}{2}
   +
   \cdots
   +
   \vat{q}{\nz - 1}
   +
   \vat{q}{\nz}

