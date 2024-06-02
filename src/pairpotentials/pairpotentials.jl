

module PairPotentials

using DiffResults, ForwardDiff, StaticArrays, Unitful
using LinearAlgebra: norm
using AtomsBase: AbstractSystem 
import AtomsCalculators

import AtomsCalculatorsUtilities.SitePotentials: SitePotential, 
                     eval_site, eval_grad_site, cutoff_radius, 
                     energy_unit, length_unit


# NOTE: this could be taken a subtype of SitePotential, but not clear 
#       that this is the best way to do it, one could gain a factor 2 
#       performance by keeping it a separate implementation. 

"""
`PairPotential`:abstractsupertype for pair potentials
"""
abstract type PairPotential <: SitePotential end


struct SimplePairPotential{ID, TC, TE} <: PairPotential where {TC<:Unitful.LengthUnits, TE<:Unitful.EnergyUnits }
   f::Function
   atom_ids::NTuple{2,ID}
   cutoff::TC
   "zero used to define output type and unit"
   zero_energy::TE  
end


struct ParametricPairPotential{ID, TP, TC, TE} <: PairPotential where {TP<:AbstractArray, TC<:Unitful.LengthUnits, TE<:Unitful.EnergyUnits }
   f::Function
   parameters::TP
   atom_ids::NTuple{2,ID}
   cutoff::TC
   "zero used to define output type and unit"
   zero_energy::TE  
end


# -------------------------------------------------------
# generic interface functionality 

# NOTE: all this needs to be cleaned up once the new 
#       AtomsCalculators interface is in place. 


cutoff_radius(pp::PairPotential) = pp.cutoff

energy_unit(pp::PairPotential) = unit(pp.zero_energy)

length_unit(pp::PairPotential) = unit(pp.cutoff)


##



function eval_site(spp::SimplePairPotential, Rs, Zs, z0)
    TR = eltype(eltype(Rs))
    if ! (z0 in spp.atom_ids)
        return zero(TR) 
    end
    tmp = zero(TR)
    id = z0 == spp.atom_ids[1] ?  spp.atom_ids[2] : spp.atom_ids[1]
    for (i, R) in zip(Zs, Rs)
        if i == id
            r = norm(R)
            tmp += spp.f(r)
        end
    end
    return tmp/2 # divide by to to get correct double count
end


function eval_grad_site(spp::SimplePairPotential, Rs, Zs, z0)
   TR = eltype(eltype(Rs))
    @assert length(Rs) == length(Zs)
    f = zeros(SVector{3, TR}, length(Zs))
    e_tmp = zero(TR)
    if ! (z0 in spp.atom_ids)  # potential is not defined for this case
        return e_tmp, f  # return zeros - this is not the optimal but will do for now
    end
    id = z0 == spp.atom_ids[1] ?  spp.atom_ids[2] : spp.atom_ids[1]
    d_result = DiffResults.DiffResult(e_tmp, e_tmp)
    for (i, Z, R) in zip(1:length(Zs), Zs, Rs)
        if Z == id
            r = norm(R)
            d_result = ForwardDiff.derivative!(d_result, spp.f, r)
            te::Float64 = DiffResults.value(d_result)  # type instablity here
            e_tmp += te
            tmp::Float64 =  DiffResults.derivative(d_result)  # type instablity here
            f[i] = ( tmp / (2r) ) * R  # divide with two here to take off double count
        end
    end
    return e_tmp/2, f
end

##



function eval_site(ppp::ParametricPairPotential, Rs, Zs, z0)
   TR = eltype(eltype(Rs))
    if ! (z0 in ppp.atom_ids)
        return zero(TR)
    end
    tmp = zero(TR)
    id = z0 == ppp.atom_ids[1] ?  ppp.atom_ids[2] : ppp.atom_ids[1]
    for (i, R) in zip(Zs, Rs)
        if i == id
            r = norm(R)
            tmp += ppp.f(r, ppp.parameters)
        end
    end
    return tmp/2 # divide by to to get correct double count
end


function eval_grad_site(ppp::ParametricPairPotential, Rs, Zs, z0)
    @assert length(Rs) == length(Zs)
    T = eltype(eltype(Rs))
    f = zeros(SVector{3, T}, length(Zs))
    e_tmp = ustrip(ppp.zero_energy)
    if ! (z0 in ppp.atom_ids)  # potential is not defined for this case
        return e_tmp, f  # return zeros - this is not the optimal but will do for now
    end
    id = z0 == ppp.atom_ids[1] ?  ppp.atom_ids[2] : ppp.atom_ids[1]
    d_result = DiffResults.DiffResult(zero(T), zero(T))
    for (i, Z, R) in zip(1:length(Zs), Zs, Rs)
        if Z == id
            r = norm(R)
            d_result = ForwardDiff.derivative!(d_result, x->ppp.f(x, ppp.parameters), r)
            te::T= DiffResults.value(d_result)  # type instablity here
            e_tmp += te
            tmp::T =  DiffResults.derivative(d_result)  # type instablity here
            f[i] = ( tmp / (2r) ) * R  # divide with two here to take off double count
        end
    end
    return e_tmp/2, f
end


# Parameter estimation ∂f/∂params
function eval_site(ppp::ParametricPairPotential, params::AbstractArray, Rs, Zs, z0)
    T = promote_type( (typeof ∘ ustrip ∘ zero)(ppp) , (eltype ∘ eltype)(Rs) )
    tmp = zeros( T, length(params) )
    if ! (z0 in ppp.atom_ids)
        return tmp
    end
    id = z0 == ppp.atom_ids[1] ?  ppp.atom_ids[2] : ppp.atom_ids[1]
    for (i, R) in zip(Zs, Rs)
        if i == id
            r = norm(R)
            tmp += ForwardDiff.gradient( a -> ppp.f(r, a), params)
        end
    end
    return tmp/2 # divide by to to get correct double count
end


function eval_grad_site(ppp::ParametricPairPotential, params::AbstractArray, Rs, Zs, z0)
    @assert length(Rs) == length(Zs)
    # Plan is to calculate ∂E/∂rᵢ∂pⱼ with Hessian calculation
    # using [r, params...] as imput.
    T = promote_type( (typeof ∘ ustrip ∘ zero)(ppp) , (eltype ∘ eltype)(Rs) )
    m = SMatrix{length(Rs[1]), length(params)}(zeros(T, length(Rs[1]), length(params)))
    f = fill( m, length(Zs) )
    if ! (z0 in ppp.atom_ids)  # potential is not defined for this case
        return f  # return zeros - this is not the optimal but will do for now
    end
    id = z0 == ppp.atom_ids[1] ?  ppp.atom_ids[2] : ppp.atom_ids[1]

    for (i, Z, R) in zip(1:length(Zs), Zs, Rs)
        if Z == id
            r = norm(R)
            hess = ForwardDiff.hessian( a -> ppp.f(a[1], a[2:end]), [r, params...] )
            f[i] = [ ( tmp / (2r) ) * Rᵢ for Rᵢ in R, tmp in @view hess[2:end, 1] ]
        end
    end
    return f
end


end # end module PairPotentials