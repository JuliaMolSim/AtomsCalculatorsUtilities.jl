# AtomsCalculatorsUtilities

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamolsim.github.io/AtomsCalculatorsUtilities.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamolsim.github.io/AtomsCalculatorsUtilities.jl/dev/)
[![Build Status](https://github.com/juliamolsim/AtomsCalculatorsUtilities.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/juliamolsim/AtomsCalculatorsUtilities.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides a variety of utilities to implement interatomic potentials, including e.g. generic functionality for site potentials and pair potentials, combination calculators, subsystem calculators and so forth. It also provides some simple tools for testing.

At the moment this package only implements the high-level AtomsCalculators 0.2 interface, but not yet the low-level interface, in particular the conventions for parameterized models still needs to be implemented. 

## List of current features

- Pair- and site potential building tools
- Utility calculators
- [i-PI](https://ipi-code.org) interface
- Testing tools for calculators

## i-PI interface

The interface is based on server and calculator, called driver, interaction.
Server is setup by either network mode or unixsocket mode.

Julia server can be started with

```julia
using AtomsCalculatorsUtilities.IPI

# network mode on localhost 
calc = IPIcalculator(ip"127.0.0.1"; port=33415)

# unixsocket mode with socket at /tmp/ipi_mysocket
calc = IPIcalculator("mysocket"; unixsocket=true,  basename="/tmp/ipi_")
```

Julia driver can be started with

```julia
using AtomsCalculatorsUtilities.IPI

sys = # generate system that sets up atom types for the calculator
calc = # generate a AtomsCalculators calculator

# network mode on localhost 
run_driver(ip"127.0.0.1", calc, sys; port=33415)

# unixsocket mode with socket at /tmp/ipi_mysocket
run_driver("mysocket", calc, sys; unixsocket=true,  basename="/tmp/ipi_")
```

Note, that you need to give driver AtomsBase structure that defines atom types.
If you want to change either number or types of atoms, you need to create new driver.

i-PI protocol requires values for energy, forces and virial.
If you don't need virial or if your calculator does not support virial calculation,
you can wrap your calculator to a utility calculator that returns zero virials

```julia
using AtomsCalculatorsUtilities.Calculators

my_zero_vir_calc = ZeroVirialCalculator(mycalc)
```
