using AtomsCalculatorsUtilities.Calculators
using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using StaticArrays
using Test
using Unitful

@testset "AtomsCalculatorsUtilities.jl" begin
    @testset "Calculators" begin include("test_calculators.jl") end
    @testset "PairPotentials" begin include("test_pairpotential.jl") end
    # @testset "FD Tests" begin include("test_fdtests.jl") end
end
