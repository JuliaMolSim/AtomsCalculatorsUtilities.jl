
using AtomsCalculators, Unitful, AtomsBase, StaticArrays, Test 
using AtomsCalculatorsUtilities

ACT = AtomsCalculatorsUtilities.Testing

##

module DemoPairCalc
   using AtomsCalculators, AtomsBase, ForwardDiff, Unitful, StaticArrays
   using AtomsCalculators: @generate_interface 
   using LinearAlgebra: norm 
   import AtomsCalculators: energy_forces_virial, 
                            potential_energy, forces, virial 

   _ustripvecvec(A::AbstractVector{<: AbstractArray}) = [ustrip.(a) for a in A] 

   abstract type AbstractPot end 

   struct Pot <: AbstractPot end 
   struct PotFerr <: AbstractPot end 
   struct PotVerr <: AbstractPot end 



   uE = u"eV" 
   uL = u"Å"
   uL_sys = u"Å"


   _v(r) = exp( - sum(abs2, r) )
   _dv(r) = ForwardDiff.gradient(_v, r)

   function _energy(X) 
      return sum(_v(X[j] - X[i]) 
                 for i = 1:length(X), j = 1:length(X) 
                  if i != j)
   end 

   function _forces(X) 
      f = zeros(SVector{3, Float64}, length(X))
      for i = 1:length(X), j = 1:length(X)
         if i != j 
            f[i] += _dv(X[j] - X[i])
            f[j] -= _dv(X[j] - X[i])
         end
      end
      return f
   end

   function _virial(X) 
      vir = @SMatrix zeros(3,3)
      for i = 1:length(X), j = 1:length(X)
         if i != j 
            𝐫 = X[j] - X[i]
            vir -= _dv(𝐫) * 𝐫'
         end
      end
      return vir 
   end
   

   # @generate_interface  ... not working as expected 
   potential_energy(sys, calc::AbstractPot; kwargs...) = 
        _energy(_ustripvecvec(position(sys))) * uE 

   forces(sys, calc::AbstractPot; kwargs...) = 
         _forces(_ustripvecvec(position(sys))) * uE / uL

   virial(sys, calc::AbstractPot; kwargs...) = 
         _virial(_ustripvecvec(position(sys))) * uE

   forces(sys, calc::PotFerr; kwargs...) = 
         0.9 * _forces(_ustripvecvec(position(sys))) * uE / uL

   virial(sys, calc::PotVerr; kwargs...) = 
         0.9 * _virial(_ustripvecvec(position(sys))) * uE


   function random_system(Nat)
      bb = [ SA[1.0,0.0,0.0] + 0.1 * rand(SVector{3, Float64}),
             SA[0.0,1.0,0.0] + 0.1 * rand(SVector{3, Float64}),
             SA[0.0,0.0,1.0] + 0.1 * rand(SVector{3, Float64}), ] * uL_sys
      X = [ Atom(1, rand(SVector{3, Float64})*uL_sys, missing) for _ = 1:5 ]
      periodic_system(X, bb)
   end

end

D = DemoPairCalc


##
# rattle = 0.1u"Å"
for rattle in (false, 0.1u"Å")
   Nat = rand(4:8) 
   sys = D.random_system(Nat)
   calc = D.Pot()
   calcFerr = D.PotFerr()
   calcVerr = D.PotVerr()

   result = ACT.fdtest(sys, calc; rattle=rattle)
   @test result.f_result
   @test result.v_result

   result = ACT.fdtest(sys, calcFerr; rattle=rattle)
   @test !result.f_result
   @test result.v_result

   result = ACT.fdtest(sys, calcVerr; rattle=rattle)
   @test result.f_result
   @test !result.v_result
end


##

Nat = 8
sys = D.random_system(Nat)
calc = D.Pot()

# the next test should fail because the tolerance is too stringent 
result = ACT.fdtest(sys, calc; rattle = 0.1u"Å", tol = 1e-10)
@test !result.f_result
@test !result.v_result

# another test that checks what happens if we change unit 
result = ACT.fdtest(sys, calc; rattle = 0.1u"Å", h0 = 1e-12 * u"m")
@test result.f_result
@test result.v_result

# and finally a test with inconsistent units 
D.uL = u"nm"
sys = D.random_system(Nat)
h0 =  1e-13 * u"m"
@show unit(position(sys, 1)[1])
@show unit(D.forces(sys, calc)[1][1])
@show unit(h0)
result = ACT.fdtest(sys, calc; rattle = 0.1u"Å", h0 = h0)

##
