
export ZeroVirialCalculator


@doc """
    ZeroVirialCalculator{T,VU}

Calculator that returns zeros for virial.
This calculator helps to add virial calculation to a calculator
that does not have it supported. This allows to embed the calculator
to calculations that need some value for virial returned.

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
mutable struct ZeroVirialCalculator{T,VU}
    calculator::T
    virial_unit::VU
    function ZeroVirialCalculator(calc; virial_unit::Unitful.EnergyUnits=u"eV")
        @warn "Setting virial to zeros leads to errors! Use on your own risk."
        new{typeof(calc), typeof(virial_unit)}(calc, virial_unit)
    end
end


AtomsCalculators.zero_forces(sys, calc::ZeroVirialCalculator) = AtomsCalculators.zero_forces(sys, calc.calculator)
AtomsCalculators.promote_force_type(sys, calc::ZeroVirialCalculator) = AtomsCalculators.promote_force_type(sys, calc.calculator)


function AtomsCalculators.potential_energy(sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.potential_energy(sys, calc.calculator; kwargs...)
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, calc::ZeroVirialCalculator; kwargs...)
    # Ideally this would be same as energy unit, but it is not possible to access it
    return zeros(3,3) * u"eV"
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

function AtomsCalculators.calculate(f::AtomsCalculators.Forces, sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.calculate(f, sys, calc.calculator; kwargs...)
end

function AtomsCalculators.calculate(e::AtomsCalculators.Energy, sys, calc::ZeroVirialCalculator; kwargs...)
    return AtomsCalculators.calculate(e, sys, calc.calculator; kwargs...)
end