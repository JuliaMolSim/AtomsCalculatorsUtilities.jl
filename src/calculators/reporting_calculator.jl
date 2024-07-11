
export ReportingCalculator

@doc """
    generate_message(sys, calculator, calc_result; kwargs...) = calc_result

This is the default function that is called when `ReportingCalculator` collects
a message. Extending this allows you to control what is reported.

This function is ment to allow setting of global stetting. If you want to
set reporting function for an individual case, give `ReportingCalculator` keyword
`message_function=my_report` where `my_report` is function that returns your message.

If function returns `nothing` the message is ignored. You can use this to control
when message is sent. 
"""
generate_message(sys, calculator, calc_result; kwargs...) = calc_result


@doc """
    ReportingCalculator{T, TC, TF}

`ReportingCalculator` collects information during calculation
and sent it to a `Channel` that can be read.

# Fields

- `calculator::T`          : caculator used in calculations
- `channel::Channel{TC}`   : `Channel` where message is put
- `message::TF`            : function that generates the message

# Creation

```julia
rcalc = ReportingCalculator(calculator, Channel(32))
rcalc = ReportingCalculator(calculator, Channel(32); message_function=my_message_function)
```

When `message_function` is omitted, `generate_message` function is used. See it for more details on how to control generated messages.

You can access the channel by calling calculator directly with `fetch` or `take!`.
"""
mutable struct ReportingCalculator{T}
    calculator::T
    channel::AbstractChannel
    message::Function
    function ReportingCalculator(
        calc, 
        channel::AbstractChannel=Channel(); 
        message_function=nothing
    )
        message = something(message_function, generate_message)
        new{typeof(calc)}(calc, channel, message)
    end
end


function Base.show(io::IO, ::MIME"text/plain", calc::ReportingCalculator)
    print(io, "ReportingCalculator")
end

Base.fetch(rcalc::ReportingCalculator) = fetch(rcalc.channel)
Base.take!(rcalc::ReportingCalculator) = take!(rcalc.channel)

#AtomsCalculators.zero_forces(sys, calc::ReportingCalculator) = AtomsCalculators.zero_forces(sys, calc.calculator)
AtomsCalculators.promote_force_type(sys::AtomsBase.AbstractSystem, calc::ReportingCalculator) = AtomsCalculators.promote_force_type(sys, calc.calculator)

AtomsCalculators.energy_unit(calc::ReportingCalculator) = AtomsCalculators.energy_unit(calc.calculator)
AtomsCalculators.length_unit(calc::ReportingCalculator) = AtomsCalculators.length_unit(calc.calculator)


function AtomsCalculators.potential_energy(
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    e = AtomsCalculators.potential_energy(sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, e; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return e
end


function AtomsCalculators.virial(
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    v = AtomsCalculators.virial(sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, v; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return v
end


function AtomsCalculators.forces(
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    f = AtomsCalculators.forces(sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, f; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return f
end


function AtomsCalculators.forces!(
    f,
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    fout = AtomsCalculators.forces!(f, sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, fout; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return fout
end


function AtomsCalculators.calculate(
    calc_method::Union{
        AtomsCalculators.Energy,
        AtomsCalculators.Forces,
        AtomsCalculators.Virial
    },
    sys, 
    calc::ReportingCalculator,
    pr=nothing,
    st=nothing; 
    kwargs...
)
    tmp = AtomsCalculators.calculate(calc_method, sys, calc.calculator, pr, st; kwargs...)
    mess = calc.message(sys, calc.calculator, tmp; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return tmp
end



function AtomsCalculators.energy_forces(
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    ef = AtomsCalculators.energy_forces(sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, ef; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return ef
end


function AtomsCalculators.energy_forces_virial(
    sys, 
    calc::ReportingCalculator; 
    kwargs...
)
    efv = AtomsCalculators.energy_forces_virial(sys, calc.calculator; kwargs...)
    mess = calc.message(sys, calc.calculator, efv; kwargs...)
    if ! isnothing(mess)
        put!(calc.channel, mess)
    end
    return efv
end