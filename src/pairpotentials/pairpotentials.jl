

module PairPotentials

using ForwardDiff, StaticArrays, Unitful
using LinearAlgebra: norm
using AtomsBase: AbstractSystem 
import AtomsCalculators
import ForwardDiff: Dual

import AtomsCalculatorsUtilities.SitePotentials
import AtomsCalculatorsUtilities.SitePotentials: SitePotential, 
                     eval_site, eval_grad_site
                     #  , cutoff_radius, energy_unit, length_unit

export SimplePairPotential


# NOTE: this could be taken a subtype of SitePotential, but not clear 
#       that this is the best way to do it, one could gain a factor 2 
#       performance by keeping it a separate implementation. 

"""
`PairPotential`:abstractsupertype for pair potentials. A concrete 
type must implement the following methods: 
- `eval_pair(V, r, [z1, z2, [ps, st]]) -> val`
- `cutoff_radius(V) -> rcut`
- `energy_unit(V) -> uE`
- `length_unit(V) -> uL`

There is a default implementation of `eval_grad_pair` that uses 
ForwardDiff. Optionally this can be overloaded e.g. if the implementation 
of `eval_pair` does not allow for Dual numbers. 

AtomsCalculators uses the Lux interface is used manage parameterized potentials. 
"""
abstract type PairPotential <: SitePotential end

""" 
    eval_pair(V, r, [z1, z2, [ps, st]]) -> val 

Evaluate a pair potential `V` acting between two particles of 
species `z1` and `z2` at a distance `r`, with parameters `ps` and state `st`.
""" 
function eval_pair end 


""" 
    eval_grad_pair(V, r, [z1, z2, [ps, st]]) -> (val, deriv)

Evaluate a pair potential `V` and it's derivative acting between two particles of 
species `z1` and `z2` at a distance `r`, with parameters `ps` and state `st`.
""" 
function eval_grad_pair(V, r, args...) 
    dr = Dual(r, 1)
    dV = eval_pair(V, dr, args...)
    return dV.value, dV.partials.values[1]
end


# ----------- Default implementation of a pair potential 
#             as a simple site potential 
# In the future with better neighbourlists this can surely 
# be significantly improved upon. 

function eval_site(V::PairPotential, Rs, Zs, z0, args...)

    @assert length(Rs) == length(Zs)

    if length(Rs) > 0 
        return sum( eval_pair(V, norm(𝐫), zj, z0, args...)
                    for (𝐫, zj) in zip(Rs, Zs) ) / 2
    else 
        # the floating point type we assume for the potential
        # if we can't infer it from eval_pair since Rs is empty...
        TR = eltype(eltype(Rs))   
        return zero(TR)
    end 

end


function eval_grad_site(V::PairPotential, Rs, Zs, z0, args...)

    @assert length(Rs) == length(Zs)
    TR = eltype(eltype(Rs)) 

    # if the input is empty, return zeros and empty forces 
    if length(Rs) == 0 
        TF = eltype(Rs) 
        return zero(TR), TF[] 
    end 

    # We can infer the actual types by evaluating the potential.
    r1 = norm(Rs[1]) 
    v, dv = eval_grad_pair(V, r1, Zs[1], z0, args...)
    TE = typeof(v) 
    TF = promote_type(TR, typeof(dv)) 

    # allocate memory for the site gradient 
    Ei = v / 2
    ∇Ei = zeros(SVector{3, TF}, length(Rs))
    ∇Ei[1] = (dv/2) * Rs[1] / r1

    for j = 2:length(Rs) 
        rj = norm(Rs[j])
        v, dv = eval_grad_pair(V, rj, Zs[j], z0, args...)
        Ei += v/2
        ∇Ei[j] = (dv/2) * Rs[j] / rj
    end

    return Ei, ∇Ei
end



## Simple Pairpotential implementation


"""
    SimplePairPotential{TF, TS, UE, UL} <: PairPotential

Create simple analytical pair potential.

To create a new potential you need to supply

- scalar to scalar function for pair energy
- two identities for atoms to for a pair (currently atomic numbers)
- cutoff radius

You can also give keywords for energy and length units that the resulting
energy, forces and virial have. Defaults to eV and Å.

# Fields
- `f::TF`                      : function that takes distance and returns energy
- `atom_species::Tuple{TS,TS}` : pair of atom species among which the potential has valid
- `rcut::Float64`              : cutoff radius that has unit `UL` 

# Example
```julia
V = SimplePairPotential(
    r -> -1/r^6
    1,          # hydrogen
    2,          # helium
    6.0u"Å";
    energy_unit=u"eV",
    length_unit=u"Å"
)
```
"""
struct SimplePairPotential{TF, TS, UE, UL} <: PairPotential
    f::TF
    atom_species::Tuple{TS,TS}
    rcut::Float64
    function SimplePairPotential(
        f, 
        z1::TS, 
        z2::TS, 
        rcut::Unitful.Length; 
        energy_unit::Unitful.EnergyUnits=u"eV", 
        length_unit::Unitful.LengthUnits=u"Å"
    )   where {TS}
        new{typeof(f), TS, energy_unit, length_unit}(f, (z1,z2), ustrip(length_unit,rcut))
    end
end

AtomsCalculators.energy_unit(::SimplePairPotential{TF,TI,UE,UL}) where{TF,TI,UE,UL} = UE
AtomsCalculators.length_unit(::SimplePairPotential{TF,TI,UE,UL}) where{TF,TI,UE,UL} = UL
SitePotentials.cutoff_radius(V::SimplePairPotential) = V.rcut * AtomsCalculators.length_unit(V)


function PairPotentials.eval_pair(V::SimplePairPotential, r, z1, z2)
    if V.atom_species == (z1, z2) || V.atom_species == (z2, z1)
        return V.f(r)
    else
        return zero(r)
    end
end


end # end module PairPotentials