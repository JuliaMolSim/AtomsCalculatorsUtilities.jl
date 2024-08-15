
export SubSystemCalculator

@doc """
    SubSystemCalculator{T, TC}

Submits subsystem to given calculator.

The purpose of this calculator is that you can split system to smaller
system that can then be calculated with e.g. with different methods.
One possible use case here is QM/MM calculations where you can split
QM system out.

The structrure is mutable to allow mutable calculators.

# Fields
- `calculator::T`  :  calculator which is used for the subsystem calculation
- `subsys::TC`     :  definition of subsystem like array of indices - has to be iterable
"""
mutable struct SubSystemCalculator{T, TC} # Mutable struct so that calculator can mutate inself
    calculator::T
    subsys::TC
    function SubSystemCalculator(calc, subsys)
        @assert applicable(first, subsys) "subsys is not iterable"
        new{typeof(calc), typeof(subsys)}(calc, subsys)
    end
end

function Base.show(io::IO, ::MIME"text/plain", calc::SubSystemCalculator)
    print(io, "SubSystemCalculator - subsystem size = ", length(calc.subsys))
end


AtomsCalculators.zero_forces(sys, calc::SubSystemCalculator) = AtomsCalculators.zero_forces(sys, calc.calculator)
AtomsCalculators.promote_force_type(sys::AtomsBase.AbstractSystem, calc::SubSystemCalculator) = AtomsCalculators.promote_force_type(sys, calc.calculator)

AtomsCalculators.energy_unit(calc::SubSystemCalculator) = AtomsCalculators.energy_unit(calc.calculator)
AtomsCalculators.length_unit(calc::SubSystemCalculator) = AtomsCalculators.length_unit(calc.calculator)

function _generate_subsys(sys, calc::SubSystemCalculator)
    @assert length(sys) >= length(calc.subsys)
    sub_atoms = [ sys[i] for i in calc.subsys  ]
    sub_sys = FlexibleSystem(
        sub_atoms;
        [ k => sys[k] for k in keys(sys) ]...
    )
    return sub_sys
end


AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(sys, calc::SubSystemCalculator; kwargs...)
    sub_sys = _generate_subsys(sys, calc)
    return AtomsCalculators.potential_energy(sub_sys, calc.calculator; kwargs...)
end


AtomsCalculators.@generate_interface function AtomsCalculators.forces!(f, sys, calc::SubSystemCalculator; kwargs...)
    @assert length(f) == length(sys)
    sub_sys = _generate_subsys(sys, calc)
    tmp_f = AtomsCalculators.zero_forces(sub_sys, calc)
    AtomsCalculators.forces!(tmp_f, sub_sys, calc.calculator; kwargs...)
    #TODO this wont work for GPU Arrays
    for (i, val) in zip(calc.subsys, tmp_f)
        f[i] += val
    end
    return f
end

AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, calc::SubSystemCalculator; kwargs...)
    sub_sys = _generate_subsys(sys, calc)
    return AtomsCalculators.virial(sub_sys, calc.calculator; kwargs...)
end


function AtomsCalculators.energy_forces(sys, calc::SubSystemCalculator; kwargs...)
    sub_sys = _generate_subsys(sys, calc)
    tmp_ef = AtomsCalculators.energy_forces(sub_sys, calc.calculator; kwargs...)
    tmp_f = tmp_ef.forces
    f = zeros(AtomsCalculators.promote_force_type(sys, calc), length(sys))
    for (i, val) in zip(calc.subsys, tmp_f)
        f[i] += val
    end
    return (energy=tmp_ef.energy, forces=f)
end

function AtomsCalculators.energy_forces_virial(sys, calc::SubSystemCalculator; kwargs...)
    sub_sys = _generate_subsys(sys, calc)
    tmp_efv = AtomsCalculators.energy_forces_virial(sub_sys, calc.calculator; kwargs...)
    tmp_f = tmp_efv.forces
    f = zeros(AtomsCalculators.promote_force_type(sys, calc), length(sys))
    for (i, val) in zip(calc.subsys, tmp_f)
        f[i] += val
    end
    return (energy=tmp_efv.energy, forces=f, virial=tmp_efv.virial)
end

## Low-level interface specials

AtomsCalculators.get_state(scalc::SubSystemCalculator) = AtomsCalculators.get_state(scalc.calculator)
AtomsCalculators.get_parameters(scalc::SubSystemCalculator) = AtomsCalculators.get_parameters(scalc.calculator)

function AtomsCalculators.set_state!(scalc::SubSystemCalculator, st)
    tmp = AtomsCalculators.set_state!(scalc.calculator, st)
    return SubSystemCalculator(tmp, scalc.subsys)
end

function AtomsCalculators.set_parameters!(scalc::SubSystemCalculator, ps) 
    tmp = AtomsCalculators.set_parameters!(scalc.calculator, ps)
    return SubSystemCalculator(tmp, scalc.subsys)
end