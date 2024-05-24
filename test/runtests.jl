using AtomsCalculatorsUtilities.Calculators
using AtomsBase
using AtomsCalculators
using AtomsCalculators.AtomsCalculatorsTesting
using StaticArrays
using Test
using Unitful

@testset "AtomsCalculatorsUtilities.jl" begin
    # Define simple calculator for the tests

    struct MyType
    end

    AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(system, calculator::MyType; kwargs...)
        return 1.0u"eV" * length(system)
    end
    
    AtomsCalculators.@generate_interface function AtomsCalculators.virial(system, calculator::MyType; kwargs...)
        return ones(3,3) * u"eV" * length(system)
    end
    
    AtomsCalculators.@generate_interface function AtomsCalculators.forces(system, calculator::MyType; kwargs...)
        f0 = SVector(1.0u"eV/Å", 1.0u"eV/Å", 1.0u"eV/Å")
        return fill(f0, length(system))
    end

    hydrogen = isolated_system([
        :H => [0, 0, 0.]u"Å",
        :H => [0, 0, 1.]u"Å",
        :H => [4., 0, 0.]u"Å",
        :H => [4., 1., 0.]u"Å"
    ])

    @testset "SubSystemCalculator" begin
        
        sub_calc = SubSystemCalculator(MyType(), 1:2)

        test_energy_forces_virial(hydrogen, sub_calc)

        f = AtomsCalculators.zero_forces(hydrogen, sub_calc)
        f_zero = f[1]
        f_one = (ones ∘ typeof)( ustrip.(f_zero) ) * unit(f_zero[1])
        

        @test AtomsCalculators.potential_energy(hydrogen, sub_calc) == 2.0u"eV"

        AtomsCalculators.forces!(f, hydrogen, sub_calc)
        @test f[1] == f_one
        @test f[2] == f_one
        @test f[3] == f_zero
        @test f[4] == f_zero

        v = AtomsCalculators.virial(hydrogen, sub_calc)
        @test v[1,1] == 2.0u"eV"
    end

    @testset "CombinationCalculator" begin
        
        co_calc = CombinationCalculator(MyType(), MyType())

        test_energy_forces_virial(hydrogen, co_calc)
        
        e = AtomsCalculators.potential_energy(hydrogen, co_calc)
        f = AtomsCalculators.forces(hydrogen, co_calc)
        v = AtomsCalculators.virial(hydrogen, co_calc)
        e_ref = 2* AtomsCalculators.potential_energy(hydrogen, MyType())
        f_ref = 2* AtomsCalculators.forces(hydrogen, MyType())
        v_ref = 2* AtomsCalculators.virial(hydrogen, MyType())
        @test e ≈ e_ref
        @test all( f .≈ f_ref )
        @test all( v .≈ v_ref )
    end

    @testset "ReportingCalculator" begin
        rcalc = ReportingCalculator(MyType(), Channel(32))
        v = AtomsCalculators.calculate(AtomsCalculators.Virial(), hydrogen, rcalc)
        @test v == fetch(rcalc)
        @test v == take!(rcalc)
        test_energy_forces_virial(hydrogen, rcalc)
    end

    @testset "ZeroVirialCalculator" begin
        zvcalc = ZeroVirialCalculator(MyType())
        test_energy_forces_virial(hydrogen, zvcalc)
    end
end
