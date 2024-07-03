using Folds

export CombinationCalculator


@doc """
    generate_keywords

This function is called when `CombinationCalculator` is used.

Default implementation will only pass keywords forward.

The call type is AtomsBase system first then all calculators and kwargs.
This will allow you to extend based on calculator type.

# Example

```julia
function AtomsCalculatorsUtilities.generate_keywords(sys, pp1::PairPotential, pp2::PairPotential; kwargs...)
    if cutoff_radius(pp1) ≈ cutoff_radius(pp2)
        nlist = PairList(sys, cutoff_radius(pp1))
        return (; :nlist => nlist, kwargs...)
    else
        return kwargs
    end
end
```

will check that PairPotentials have same cutoff radius.
Then calculates pairlist and passes it forward as a keyword. 
"""
generate_keywords(sys, calculators...; kwargs...) = kwargs


@doc """
    CombinationCalculator{N}

You can combine several calculators to one calculator with this.
Giving keyword argument `executor=SequentialEx()` toggles on multithdeaded execution
of calculators. Using `executor=DistributedEx()` executes calculators using multiprocessing.

Other use case is editing keywords that are passed on the calculators.
E.g. you can generate new keyword argument that is then passed to all calculators.
This allows you to share e.g. a pairlist between calculators.

To control what keywords are passed you need to extend `generate_keywords` function.


# Fields

- calculators::NTuple{N,Any}  : NTuple that holds calculators
- executor::Any               : Transducers executor used to execute calculation - default SequentialEx
- keywords::Function          : function used to generate keywords for calculators

# Creation

```julia
CombinationCalculator( calc1, calc2, ...; executor=SequentialEx())
```

"""
mutable struct CombinationCalculator{N, T} # Mutable struct so that calculators can mutate themself
    calculators::NTuple{N,Any}
    executor::Any
    keywords::Function
    function CombinationCalculator(calculators...; executor=SequentialEx(), keyword_generator=nothing)
        kgen = something(keyword_generator, generate_keywords)
        new{length(calculators), typeof(kgen)}(calculators, executor, kgen)
    end
end

function Base.show(io::IO, ::MIME"text/plain", calc::CombinationCalculator)
    print(io, "CombinationCalculator - ", length(calc) , " calculators")
end

Base.length(cc::CombinationCalculator) = length(cc.calculators)

Base.getindex(cc::CombinationCalculator, i) = cc.calculators[i]
Base.lastindex(cc::CombinationCalculator) = length(cc)
Base.firstindex(cc::CombinationCalculator) = 1


AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    return Folds.sum( calc.calculators, calc.executor ) do c
        AtomsCalculators.potential_energy(sys, c; new_kwargs...)
    end
end

# We don't use AtomsCalculators.@generate_interface here
# as we want special version for forces!
function AtomsCalculators.forces(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    return Folds.sum( calc.calculators, calc.executor ) do c
        AtomsCalculators.forces(sys, c; new_kwargs...)
    end
end


function AtomsCalculators.calculate( ::AtomsCalculators.Forces, sys, calc::CombinationCalculator; kwargs...)
    f = AtomsCalculators.forces(sys, calc; kwargs...)
    return (; :forces => f)
end


function AtomsCalculators.forces!(f, sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)

    # Nonallocating forces is only truly nonallocating when sequential
    foreach( calc.calculators ) do cal
        AtomsCalculators.forces!(f, sys, cal; new_kwargs...)
    end
    return f
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    return Folds.sum( calc.calculators ) do c
        AtomsCalculators.virial(sys, c; new_kwargs...)
    end
end


function AtomsCalculators.energy_forces(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    tmp = Folds.sum( calc.calculators, calc.executor ) do c
        ef = AtomsCalculators.energy_forces(sys, c; new_kwargs...)
        [ef.energy, ef.forces]
    end
    return (energy=tmp[1], forces=tmp[2])
end

function AtomsCalculators.energy_forces_virial(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    tmp = Folds.sum( calc.calculators, calc.executor ) do c
        efv = AtomsCalculators.energy_forces_virial(sys, c; new_kwargs...)
        [efv.energy, efv.forces, efv.virial]
    end
    return (energy=tmp[1], forces=tmp[2], virial=tmp[3])
end