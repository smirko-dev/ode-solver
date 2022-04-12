# ODE Solver

![languages](https://img.shields.io/badge/languages-C++%20-blue)
[![version](https://img.shields.io/badge/version-%200.1-blue)](https://github.com/smirko-dev/ode-solver/blob/main/CHANGELOG.md)
[![](https://img.shields.io/badge/license-MIT-blue)](https://github.com/smirko-dev/ode-solver/blob/main/LICENSE)
[![Build Actions Status](https://github.com/smirko-dev/ode-solver/workflows/Build/badge.svg)](https://github.com/smirko-dev/ode-solver/actions)

## Description

The [ODESolver](ode) library comes with two interface classes and the following implementations
- Runge Kutta
- Euler
- Mid Point
- Velocity Verlet

## How to build

```sh
git clone git@github.com:smirko-dev/ode-solver.git
cmake -B build -DCMAKE_BUILD_TYPE=Release -G "Visual Studio 16 2019"
cmake --build build --config Release
```

## Examples

### [Minimal example](ode/README.md#Example)

Example implementation of an ODEFunction using the Euler ODESolver.

### [Molecular Dynamics](moleculardynamics)

With the [Lennard Jones Potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) the movements of molecules in a system can be calculated. This ODE 2nd order can be solved by using the [Velocity Verlet](https://en.wikipedia.org/wiki/Verlet_integration) algorithm.

### [Planet Dynamics](planetdynamics)

With the [Newton's law of universal gravitation](https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation) planet movement can be calculated. This ODE 1st order can be solved by using the [Runge Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) algorithm. 
