# ODE Solver

![languages](https://img.shields.io/badge/languages-C++%20-blue)
[![version](https://img.shields.io/badge/version-%200.1-blue)](https://github.com/smirko-dev/ode-solver/blob/main/CHANGELOG.md)
[![](https://img.shields.io/badge/license-MIT-blue)](https://github.com/smirko-dev/ode-solver/blob/main/LICENSE)
[![Build Actions Status](https://github.com/smirko-dev/ode-solver/workflows/Build/badge.svg)](https://github.com/smirko-dev/ode-solver/actions)

## Description

### OdeFunction

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation.

The parameters are given in a single vector. The setParams and getParams methods implement the mapping of the parameters. 

### OdeSolver

The ODE solver provides an interface for certain implementation´s. 

## Examples

### Molecular Dynamics

With the [Lennard Jones Potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) the movements of molecules in a system can be calculated which is an ODE 2nd order by using the [Velocity Verlet](https://en.wikipedia.org/wiki/Verlet_integration) algorithm.

### Planet Dynamics

With the [Newton's law of universal gravitation](https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation) planet movement can be calculated which is an ODE 1st order by using the [Runge Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) algorithm. 
