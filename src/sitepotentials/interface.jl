"""
`SitePotential`: abstractsupertype for generic site potentials. Concrete subtypes 
should overload 
```julia 
cutoff_radius
eval_site
eval_grad_site
energy_unit 
length_unit
```
Optional methods are 
```julia 
hessian_site
block_hessian_site
precon_site
block_precon_site
```
"""
abstract type SitePotential end

"""
   cutoff_radius(V) 

For a site potential calculator `V`, this function should return the
cutoff radius used to construct the neighbourlist. 
"""
function cutoff_radius end 

"""
If `V <: SitePotential` then it should implement the method
```julia 
val = eval_site(V, Rs, Zs, z0)
```
where `Rs::AbstractVector{<: SVector{3}}` and `Zs::AbstractVector` of atom ids 
(e.g. atomic numbers), while `z0` is a single atom id. 

The output `val` should be a single number, namely the site energy.
"""
function eval_site end


"""
If `V <: SitePotential` then it should implement the method
```julia 
dv = eval_grad_site(V, Rs, Zs, z0)
```
where `Rs::AbstractVector{<: SVector{3}}` and `Zs::AbstractVector` of 
atom ids (e.g., atomic numbers), while `z0` is a single atom id. 

The output `dv` should be an `AbstractVector` containing  
`SVector{3,T}` blocks.
"""
function eval_grad_site end 


function block_precon_site end 

"""
If `V <: SitePotential` then it can implement the method
```julia 
Pblock = precon(V, Rs, Zs, z0)
```
where `Rs::AbstractVector{<: SVector{3}}` and `Zs::AbstractVector` of 
atom ids (e.g., atomic numbers), while `z0` is a single atom id. 
The output `Pblock` should be an `AbstractMatrix` containing 
`SMatrix{3,3,T}` blocks. 

Unlike `eval_site` and `eval_grad_site`, this method is optional. It 
can be used to speedup geometry optimization, sampling and related tasks. 
"""
function precon_site end 

function block_hessian_site end 

function hessian_site end 

function energy_unit end 

function length_unit end 


force_unit(V::SitePotential) = energy_unit(V) / length_unit(V)


# ---------------------------------------------------------------
# some additional utilities that could be overloaded as needed 
# but where defaults are probably fine 

# default is atomic number, but one could use atomic symbol or something else
# ideally the default should be changed to ChemicalSymbol  i.e. 
# an isbits type that is not a Number. 
@inline get_id(at, V, i) = atomic_number(at, i)

function length_types(sys::AbstractSystem)
   r = position(sys, 1)[1][1]
   return eltype(ustrip(r)), unit(r)
end

function init_energy(sys::AbstractSystem, V::SitePotential) 
   TL, uL = length_types(sys)
   uE = energy_unit(V)
   return zero(TL) * uE 
end

function init_forces(sys::AbstractSystem{D}, V::SitePotential) where {D} 
   TL, uL = length_types(sys)
   return zeros(SVector{D, TL}, length(sys)) * energy_unit(V) / uL 
end 

function init_virial(sys::AbstractSystem{D}, V::SitePotential) where {D}
   TL, uL = length_types(sys)
   return zero(SMatrix{D,D,TL}) * energy_unit(V)
end
