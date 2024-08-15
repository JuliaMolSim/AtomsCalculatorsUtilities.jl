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
mutable struct CombinationCalculator{N,T,TE,TL} # Mutable struct so that calculators can mutate themself
    calculators::NTuple{N,Any}
    executor::Any
    keywords::T
    energy_unit::TE
    length_unit::TL
    function CombinationCalculator(
        calculators...; 
        executor=SequentialEx(), 
        keyword_generator=nothing,
        energy_unit=AtomsCalculators.energy_unit(calculators[1]),
        length_unit=AtomsCalculators.length_unit(calculators[1])
    )
        kgen = something(keyword_generator, generate_keywords)
        new{length(calculators), typeof(kgen), typeof(energy_unit), typeof(length_unit)}(
            calculators, 
            executor, 
            kgen, 
            energy_unit, 
            length_unit
        )
    end
end


function Base.show(io::IO, ::MIME"text/plain", calc::CombinationCalculator)
    print(io, "CombinationCalculator - ", length(calc) , " calculators")
end

Base.length(cc::CombinationCalculator) = length(cc.calculators)

Base.getindex(cc::CombinationCalculator, i) = cc.calculators[i]
Base.lastindex(cc::CombinationCalculator) = length(cc)
Base.firstindex(cc::CombinationCalculator) = 1

AtomsCalculators.energy_unit(cc::CombinationCalculator) = cc.energy_unit
AtomsCalculators.length_unit(cc::CombinationCalculator) = cc.length_unit


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    et::AtomsCalculators.Energy, 
    sys, 
    calc::CombinationCalculator, 
    pr=nothing, 
    st=nothing; 
    kwargs...
)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    if isnothing(pr)
        tpr = Tuple( nothing for _ in 1:length(calc))
    else
        tpr = Tuple( pr[i] for i in 1:length(calc)  ) # This checks for correct length too
    end
    if isnothing(st)
        tst = Tuple( nothing for _ in 1:length(calc))
    else
        tst = Tuple( st[i] for i in 1:length(calc)  ) # This checks for correct length too
    end
    tmp =  Folds.map( zip(calc.calculators, tpr, tst), calc.executor ) do (c, p, s)
        AtomsCalculators.calculate(et, sys, c, p, s; new_kwargs...)
    end
    Etot = sum( tmp ) do x
        x.energy
    end
    Etot = uconvert(calc.energy_unit, Etot)
    st_out = Tuple( x.state for x in tmp )
    return ( energy=Etot, state=st_out  )
end

# We don't use AtomsCalculators.@generate_interface here
# as we want special version for forces!
function AtomsCalculators.forces(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    f =  Folds.sum( calc.calculators, calc.executor ) do c
        AtomsCalculators.forces(sys, c; new_kwargs...)
    end
    fz = AtomsCalculators.zero_forces(sys, calc)
    fz .+= f
    return fz
end


function AtomsCalculators.calculate(
    ft::AtomsCalculators.Forces, 
    sys, 
    calc::CombinationCalculator, 
    pr=nothing, 
    st=nothing; 
    kwargs...
)   
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)

    # Check and prepare parameters and state
    if isnothing(pr)
        tpr = Tuple( nothing for _ in 1:length(calc))
    else
        tpr = Tuple( pr[i] for i in 1:length(calc)  ) # This checks for correct length too
    end
    if isnothing(st)
        tst = Tuple( nothing for _ in 1:length(calc))
    else
        tst = Tuple( st[i] for i in 1:length(calc)  ) # This checks for correct length too
    end
    
    tmp =  Folds.map( zip(calc.calculators, tpr, tst), calc.executor ) do (c, p, s)
        AtomsCalculators.calculate(ft, sys, c, p, s; new_kwargs...)
    end
    Ftot = sum( tmp ) do x
        x.forces
    end
    fz = AtomsCalculators.zero_forces(sys, calc)
    fz .+= Ftot
    st_out = Tuple( x.state for x in tmp )
    return (forces = fz, state = st_out)
end


function AtomsCalculators.forces!(f, sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)

    # Nonallocating forces is only truly nonallocating when sequential
    foreach( calc.calculators ) do cal
        AtomsCalculators.forces!(f, sys, cal; new_kwargs...)
    end
    return f
end


AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
    vt::AtomsCalculators.Virial,
    sys, 
    calc::CombinationCalculator,
    pr::Union{Nothing,Tuple}=nothing,
    st::Union{Nothing,Tuple}=nothing,; 
    kwargs...
)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)

    # Check and prepare parameters and state
    if isnothing(pr)
        tpr = Tuple( nothing for _ in 1:length(calc))
    else
        tpr = Tuple( pr[i] for i in 1:length(calc)  ) # This checks for correct length too
    end
    if isnothing(st)
        tst = Tuple( nothing for _ in 1:length(calc))
    else
        tst = Tuple( st[i] for i in 1:length(calc)  ) # This checks for correct length too
    end

    tmp =  Folds.map( zip(calc.calculators, tpr, tst), calc.executor ) do (c, p, s)
        AtomsCalculators.calculate(vt, sys, c, p, s; new_kwargs...)
    end
    
    # Gather output
    vir_tot = sum( tmp ) do x
        x.virial
    end
    Vtot = uconvert.(calc.energy_unit, vir_tot)
    st_out = Tuple( x.state for x in tmp )
    return ( virial=Vtot, state=st_out  )
end


function AtomsCalculators.energy_forces(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    tmp = Folds.sum( calc.calculators, calc.executor ) do c
        ef = AtomsCalculators.energy_forces(sys, c; new_kwargs...)
        [ef.energy, ef.forces]
    end
    tmp[1] = uconvert(calc.energy_unit, tmp[1])
    fz = AtomsCalculators.zero_forces(sys, calc)
    fz .+= tmp[2]
    tmp[2] = fz
    return (energy=tmp[1], forces=tmp[2])
end

function AtomsCalculators.energy_forces_virial(sys, calc::CombinationCalculator; kwargs...)
    new_kwargs = calc.keywords(sys, calc.calculators...; kwargs...)
    tmp = Folds.sum( calc.calculators, calc.executor ) do c
        efv = AtomsCalculators.energy_forces_virial(sys, c; new_kwargs...)
        [efv.energy, efv.forces, efv.virial]
    end
    tmp[1] = uconvert(calc.energy_unit, tmp[1])
    fz = AtomsCalculators.zero_forces(sys, calc)
    fz .+= tmp[2]
    tmp[2] = fz
    tmp[3] = uconvert.(calc.energy_unit, tmp[3])
    return (energy=tmp[1], forces=tmp[2], virial=tmp[3])
end




## Low-level interface specials

function AtomsCalculators.get_state(ccalc::CombinationCalculator)
    return Tuple( AtomsCalculators.get_state(calc) for calc in ccalc.calculators  )
end

function AtomsCalculators.get_parameters(ccalc::CombinationCalculator)
    return Tuple( AtomsCalculators.get_parameters(calc) for calc in ccalc.calculators  )
end

function AtomsCalculators.set_state!(ccalc::CombinationCalculator, states)
    tmp = map( zip(ccalc.calculators, states) ) do (calc, st)
        AtomsCalculators.set_state!(calc, st)
    end
    return CombinationCalculator(
        tmp...;
        executor=ccalc.executor, 
        keyword_generator=ccalc.keywords,
        energy_unit=ccalc.energy_unit,
        length_unit=ccalc.length_unit
    )
end


function AtomsCalculators.set_parameters!(ccalc::CombinationCalculator, parameters)
    tmp = map( zip(ccalc.calculators, parameters) ) do (calc, ps)
        AtomsCalculators.set_parameters!(calc, ps)
    end
    return CombinationCalculator(
        tmp...;
        executor=ccalc.executor, 
        keyword_generator=ccalc.keywords,
        energy_unit=ccalc.energy_unit,
        length_unit=ccalc.length_unit
    )
end
