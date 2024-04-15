#########################
Simple Bubbly Flow Solver
#########################

|License| |LastCommit| |CI| |DOCS|

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleBubblyFlowSolver
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleBubblyFlowSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleBubblyFlowSolver/commits/main

.. |CI| image:: https://github.com/NaokiHori/SimpleBubblyFlowSolver/actions/workflows/ci.yml/badge.svg

.. |DOCS| image:: https://github.com/NaokiHori/SimpleBubblyFlowSolver/actions/workflows/documentation.yml/badge.svg
.. _DOCS: https://naokihori.github.io/SimpleBubblyFlowSolver

.. image:: https://github.com/NaokiHori/SimpleBubblyFlowSolver/blob/main/docs/source/thumbnail.jpg
   :target: https://youtu.be/Xr14Kw61ByA
   :width: 100%

********
Overview
********

This library numerically solves the motion of two-liquid mixtures separated by free surfaces using finite-difference and volume-of-fluid methods.

Specifically, its aim is to simulate air-water flows characterised by significant contrasts in density and viscosity at high Reynolds numbers.

***********
Quick start
***********

Fetch source

.. code-block:: console

   git clone --recurse-submodules https://github.com/NaokiHori/SimpleBubblyFlowSolver
   cd SimpleBubblyFlowSolver

Initialise flow fields (needs ``Python3`` with ``NumPy``)

.. code-block:: console

   cd initial_condition
   make output
   sh main.sh
   cd ..

Build and run

.. code-block:: console

   make output
   make all
   sh exec/main.sh

This simulates the motion of a 2D rising bubble:

.. image:: https://github.com/NaokiHori/SimpleBubblyFlowSolver/blob/main/docs/source/sample.jpg

*************
Documentation
*************

The governing equations, the numerical methods employed, and the discretisations are briefly discussed `here <https://naokihori.github.io/SimpleBubblyFlowSolver/index.html>`_.

**********
3D version
**********

Checkout ``3d`` branch.
Initialise flow fields by yourself.

***************
Acknowledgement
***************

I would like to acknowledge Dr. Kevin Zhong for fruitful discussion regarding the viscosity formulation.

