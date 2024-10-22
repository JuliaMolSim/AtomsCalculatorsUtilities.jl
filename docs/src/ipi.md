# i-PI driver interface

[i-PI](https://github.com/i-pi/i-pi) is a universal force engine that defines
protocol for energy, force and virial calculation.

AtomsCalculatorUtilities has a submodule IPI that implement i-PI interface for
AtomsCalculators compatable caluculators.

## i-PI basic consepts

The interface is based on server and calculator, called driver, interaction.
Server is setup by either network mode or unixpipe mode.

Following starts [ASE](https://wiki.fysik.dtu.dk/ase/index.html) SocketIOCalculator that is a i-PI server implementation, using default port 31415

```python
from ase.calculators.socketio import SocketIOCalculator
calc = SocketIOCalculator(log="test.log")
```

You then need to connect a calculator (driver) to it. To do this you need to tell the calculator what kind of system (atom symbols) are used in the calculation, as i-PI protocol does not transfer atomic symbols.

After you have connected the calculator (driver), you can use the ASE calculator to perform calculations.

## Use Julia driver

To start Julia driver you need to first create a AtomsCalculators compatable calculator and an AtomsBase system structure that initializes atomic symbols.

Then you connect them to a driver by calling `run_driver`

```julia
using AtomsCalculatorsUtilities.IPI

mycalc = # Create your calculator
sys_init = # Create initial system

run_driver("127.0.0.1", mycalc, sys-init)
```

After this the server works as a calculator, using your Julia calculator as a backend.

## Things to note

i-PI protocol will always calculate energy, forces and virial, thus your calculator need to support virial
calculation.

For calculators that do not support virial calculations `AtomsUtilityCalcualtors.Calculators` implements
`ZeroVirialCalculator` that adds virial caluculation support for your calculator. Note that this is in correct. So use it at your own risk!

One good use case for this is constant volume simulations, where
virial is ignored, but i-PI protocol still calculates it.

## Working example

Start SocketIOCalculator in python

```python
from ase import Atoms
from ase.calculators.socketio import SocketIOCalculator

h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.7]])

calc = SocketIOCalculator(log="test.log")

h2.calc = calc
h2.get_forces()  # this will stall untill you connect a driver
```

In Julia you then create a driver

```julia
using AtomsCalculatorsUtilities.Calculators
using AtomsCalculatorsUtilities.IPI
using ORCAcalculator
using Unitful
using AtomsBase

ox = ORCAexecutable()
om = ORCAmethod("blyp def2-TZVPP TIGHTSCF defgrid3")

hydrogen = isolated_system([
    :H => [0, 0, 0.]u"Å",
    :H => [0, 0, 1.]u"Å"
])

orca = ORCAcalculatorbundle(ox, om)

calc = ZeroVirialCalculator(orca)

run_driver("127.0.0.1", calc, hydrogen)
```

Now you can perform calculations in ASE.
