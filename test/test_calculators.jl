using AtomsCalculatorsUtilities.Calculators
using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using StaticArrays
using Test
using Unitful


struct MyType
end

AtomsCalculators.energy_unit(::MyType) = u"eV"
AtomsCalculators.length_unit(::MyType) = u"Å"

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

    AtomsCalculators.Testing.test_energy_forces_virial(hydrogen, sub_calc)

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

    st = AtomsCalculators.get_state(sub_calc)
    ps = AtomsCalculators.get_parameters(sub_calc)
    nscalc = AtomsCalculators.set_state!(sub_calc, st)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nscalc)
    nscalc = AtomsCalculators.set_parameters!(sub_calc, ps)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nscalc)
end


@testset "CombinationCalculator" begin
    
    co_calc = CombinationCalculator(MyType(), MyType())

    AtomsCalculators.Testing.test_energy_forces_virial(hydrogen, co_calc)
    
    e = AtomsCalculators.potential_energy(hydrogen, co_calc)
    f = AtomsCalculators.forces(hydrogen, co_calc)
    v = AtomsCalculators.virial(hydrogen, co_calc)
    e_ref = 2* AtomsCalculators.potential_energy(hydrogen, MyType())
    f_ref = 2* AtomsCalculators.forces(hydrogen, MyType())
    v_ref = 2* AtomsCalculators.virial(hydrogen, MyType())
    @test e ≈ e_ref
    @test all( f .≈ f_ref )
    @test all( v .≈ v_ref )

    st = AtomsCalculators.get_state(co_calc)
    ps = AtomsCalculators.get_parameters(co_calc)
    nccalc = AtomsCalculators.set_state!(co_calc, st)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nccalc)
    nscalc = AtomsCalculators.set_parameters!(co_calc, ps)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nccalc)
end

@testset "ReportingCalculator" begin
    rcalc = ReportingCalculator(MyType(), Channel(32))
    v = AtomsCalculators.calculate(AtomsCalculators.Virial(), hydrogen, rcalc)
    @test v == fetch(rcalc)
    @test v == take!(rcalc)
    AtomsCalculators.Testing.test_energy_forces_virial(hydrogen, rcalc)

    st = AtomsCalculators.get_state(rcalc)
    ps = AtomsCalculators.get_parameters(rcalc)
    nrcalc = AtomsCalculators.set_state!(rcalc, st)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nrcalc)
    nrcalc = AtomsCalculators.set_parameters!(rcalc, ps)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nrcalc)
end

@testset "ZeroVirialCalculator" begin
    zvcalc = ZeroVirialCalculator(MyType())
    v = AtomsCalculators.virial(hydrogen, zvcalc)
    @test all(iszero, v)
    AtomsCalculators.Testing.test_energy_forces_virial(hydrogen, zvcalc)

    st = AtomsCalculators.get_state(zvcalc)
    ps = AtomsCalculators.get_parameters(zvcalc)
    nzcalc = AtomsCalculators.set_state!(zvcalc, st)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nzcalc)
    nzcalc = AtomsCalculators.set_parameters!(zvcalc, ps)
    AtomsCalculators.Testing.test_potential_energy(hydrogen, nzcalc)
end
