####################
Stream-wise operator
####################

*******
Average
*******

Average at stream-wise cell faces:

.. math::

   \vat{
      \ave{
         q
      }{
         \gy
      }
   }{
      j + \frac{1}{2}
   }
   =
   \frac{1}{2} \vat{q}{j} + \frac{1}{2} \vat{q}{j + 1}

Average at stream-wise cell centers:

.. math::

   \vat{
      \ave{
         q
      }{
         \gy
      }
   }{
      j
   }
   =
   \frac{1}{2} \vat{q}{j - \frac{1}{2}} + \frac{1}{2} \vat{q}{j + \frac{1}{2}}

***************
Differentiation
***************

Differentiation at stream-wise cell faces:

.. math::

   \vat{
      \dif{
         q
      }{
         \gy
      }
   }{
      j + \frac{1}{2}
   }
   =
   -
   \vat{q}{j}
   +
   \vat{q}{j + 1}

Differentiation at stream-wise cell centers:

.. math::

   \vat{
      \dif{
         q
      }{
         \gy
      }
   }{
      j
   }
   =
   -
   \vat{q}{j - \frac{1}{2}}
   +
   \vat{q}{j + \frac{1}{2}}

*********
Summation
*********

Summation of a quantity defined at stream-wise cell faces:

.. math::

   \sumyf q
   \equiv
   \vat{q}{\frac{1}{2}}
   +
   \vat{q}{\frac{3}{2}}
   +
   \cdots
   +
   \vat{q}{\ny - \frac{3}{2}}
   +
   \vat{q}{\ny - \frac{1}{2}}

Summation of a quantity defined at stream-wise cell centers:

.. math::

   \sumyc q
   \equiv
   \vat{q}{1}
   +
   \vat{q}{2}
   +
   \cdots
   +
   \vat{q}{\ny - 1}
   +
   \vat{q}{\ny}

