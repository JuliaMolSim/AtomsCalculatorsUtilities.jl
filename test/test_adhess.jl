

using AtomsBase, Unitful, ExtXYZ, AtomsCalculators, AtomsBuilder,
      EmpiricalPotentials, StaticArrays, Test, JSON, ForwardDiff 
using AtomsCalculators.AtomsCalculatorsTesting
using LinearAlgebra: dot, norm, I 
using EmpiricalPotentials: cutoff_radius, StillingerWeber
using AtomsCalculators: potential_energy, forces

EP = EmpiricalPotentials
ACT = AtomsCalculators.AtomsCalculatorsTesting

sw = StillingerWeber()

function rand_Si_struct(nrep = 2, r = 0.1) 
   rattle!(bulk(:Si, cubic=true), r)   
end

function rand_Si_env(args...)
   sys = rand_Si_struct(args...)
   nlist = EmpiricalPotentials.PairList(sys, cutoff_radius(StillingerWeber()))
   Js, Rs, Zs, z0 = EmpiricalPotentials.get_neighbours(sys, sw, nlist, 1)
   return Rs, Zs, z0 
end

##

# preliminary finit-difference test, this is to be incorporated into the 
# FD tests in the test suite in AtomsCalculators

for ntest = 1:10 
   Rs, Zs, z0 = rand_Si_env()
   Us = randn(eltype(Rs), length(Rs))
   Vs = randn(eltype(Rs), length(Rs))
   F(t) = dot( EP.eval_grad_site(sw, Rs + t * Us, Zs, z0)[2], Vs )
   dF0 = dot( EP.ad_block_hessian_site(sw, Rs, Zs, z0) * Us, Vs )

   @test ACT.fdtest(F, t -> dF0, 0.0; verbose = false )
end

##
# a quick test that shows the ForwardDiff implementation is quite decent.
# using BenchmarkTools
# Rs, Zs, z0 = rand_Si_env()
# @info("Timing of eval_grad_site")
# @btime EP.eval_grad_site($sw, $Rs, $Zs, $z0)
# @info("Timing of ad_site_hessian")
# @btime EP.ad_site_hessian($sw, $Rs, $Zs, $z0)

## 
# look at the global hessian now 
# first in block format 

for ntest = 1:10 
   sys = rand_Si_struct()
   H = EP.block_hessian(sys, sw)

   X0 = position(sys)
   uL = unit(X0[1][1])
   Us = randn(SVector{3, Float64}, length(X0)) * uL 
   Vs = randn(SVector{3, Float64}, length(X0)) * uL / u"eV" 
   _sys(X) = FastSystem(bounding_box(sys), 
                        boundary_conditions(sys), 
                        X, 
                        atomic_symbol(sys), 
                        atomic_number(sys), 
                        atomic_mass(sys))

   F(t) = - dot(forces(_sys(X0 + t * Us), sw), Vs) 
   F(0.0)
   H0 = EP.block_hessian(_sys(X0), sw)
   dF0 = dot(H0 * Us, Vs)
   @test ACT.fdtest(F, t -> dF0, 0.0; verbose = false )
end

##
# now global hessian in blas format 

function _globify(Abl::AbstractMatrix{SMatrix{D, D, T}}) where {D, T} 
   A = zeros(T, D * size(Abl, 1), D * size(Abl, 2))
   for j1 = 1:size(Abl, 1), j2 = 1:size(Abl, 2) 
      J1 = D * (j1-1) .+ (1:D)
      J2 = D * (j2-1) .+ (1:D)
      A[J1, J2] .= Abl[j1, j2]
   end
   return A 
end

for ntest = 1:10 
   Rs, Zs, z0 = rand_Si_env()
   Hi = EP.hessian_site(sw, Rs, Zs, z0)
   Hi_bl = EP.block_hessian_site(sw, Rs, Zs, z0)
   @test Hi ≈ _globify(Hi_bl)
end

for ntest = 1:10 
   sys = rand_Si_struct()
   H = EP.hessian(sys, sw)
   Hbl = EP.block_hessian(sys, sw)
   @test H ≈ _globify(Hbl)
end
