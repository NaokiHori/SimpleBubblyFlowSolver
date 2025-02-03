# Simple Bubbly Flow Solver

![License](https://img.shields.io/github/license/NaokiHori/SimpleBubblyFlowSolver) [License](https://opensource.org/licenses/MIT)
![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SimpleBubblyFlowSolver/main) [Last Commit](https://github.com/NaokiHori/SimpleBubblyFlowSolver/commits/main)
![CI](https://github.com/NaokiHori/SimpleBubblyFlowSolver/actions/workflows/ci.yml/badge.svg)
![DOCS](https://github.com/NaokiHori/SimpleBubblyFlowSolver/actions/workflows/documentation.yml/badge.svg) [Documentation](https://naokihori.github.io/SimpleBubblyFlowSolver)

[![Thumbnail](https://github.com/NaokiHori/SimpleBubblyFlowSolver/blob/main/docs/source/thumbnail.jpg)](https://youtu.be/Xr14Kw61ByA)

## Overview

This library numerically solves the motion of two-liquid mixtures separated by free surfaces using finite-difference and volume-of-fluid methods.

Specifically, it aims to simulate air-water flows characterized by contrasts in density and viscosity at high (in terms of the direct numerical simulations) Reynolds numbers.

## Quick Start

### Fetch Source

```console
git clone --recurse-submodules https://github.com/NaokiHori/SimpleBubblyFlowSolver
cd SimpleBubblyFlowSolver
```

### Initialize Flow Fields (Requires `Python3` with `NumPy`)

```console
cd initial_condition
make output
sh main.sh
cd ..
```

### Build and Run

```console
make output
make all
sh exec/main.sh
```

This simulates the motion of a 2D rising bubble:

![Sample](https://github.com/NaokiHori/SimpleBubblyFlowSolver/blob/main/docs/source/sample.jpg)

## Documentation

The governing equations, numerical methods employed, and discretizations are briefly discussed [here](https://naokihori.github.io/SimpleBubblyFlowSolver/index.html).

## 3D Version

Check out the `3d` branch. 
You will need to initialize the flow fields manually.

## Acknowledgement

I would like to acknowledge Dr. Kevin Zhong for fruitful discussions regarding the viscosity formulation.

