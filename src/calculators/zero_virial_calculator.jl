
export ZeroVirialCalculator

import AtomsCalculators: energy_unit 

@doc """
    ZeroVirialCalculator{T,VU}

Calculator that returns zeros for virial.
This calculator helps to add virial calculation to a calculator
that does not have it supported. This allows to embed the calculator
to calculations that need some value for virial returned.

Note that forcing zero virial is incorrect. Use at your own risk.

# Fields
- calculator::T    :  calculator used for potential energy and forces
- virial_unit::VU  :  unit given to zero virial - default eV

# Creation

```julia
ZeroVirialCalculator(
    mycalculator;
    virial_unit::Unitful.EnergyUnits=u"eV"
)
```
"""
mutable struct ZeroVirialCalculator{TCALC}
    calculator::TCALC
    function ZeroVirialCalculator(calc)
        @warn "Setting virial to zeros leads to errors! Use on your own risk."
        new{typeof(calc)}(calc)
    end
end


#AtomsCalculators.zero_forces(sys, calc::ZeroVirialCalculator) = AtomsCalculators.zero_forces(sys, calc.calculator)
AtomsCalculators.promote_force_type(sys::AtomsBase.AbstractSystem, calc::ZeroVirialCalculator) = AtomsCalculators.promote_force_type(sys, calc.calculator)

AtomsCalculators.energy_unit(calc::ZeroVirialCalculator) = AtomsCalculators.energy_unit(calc.calculator)
AtomsCalculators.length_unit(calc::ZeroVirialCalculator) = AtomsCalculators.length_unit(calc.calculator)


function AtomsCalculators.potential_energy(sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.potential_energy(sys, calc.calculator; kwargs...)
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, calc::ZeroVirialCalculator; kwargs...)
    D = AtomsBase.n_dimensions(sys) 
    return zeros(D, D) * energy_unit(calc)
end


function AtomsCalculators.energy_forces(sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.energy_forces(sys, calc.calculator; kwargs...)
end

function AtomsCalculators.forces(sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.forces(sys, calc.calculator; kwargs...)
end

function AtomsCalculators.forces!(f, sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.forces!(f, sys, calc.calculator; kwargs...)
end

function AtomsCalculators.calculate(
    f::AtomsCalculators.Forces, 
    sys, 
    calc::ZeroVirialCalculator,
    pr=nothing,
    st=nothing; 
    kwargs...
)
    return AtomsCalculators.calculate(f, sys, calc.calculator, pr, st; kwargs...)
end

function AtomsCalculators.calculate(
    e::AtomsCalculators.Energy, 
    sys, 
    calc::ZeroVirialCalculator,
    pr=nothing,
    st=nothing; 
    kwargs...
)
    return AtomsCalculators.calculate(e, sys, calc.calculator; kwargs...)
end