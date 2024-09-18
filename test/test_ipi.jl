using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using AtomsCalculatorsUtilities.IPI
using AtomsCalculatorsUtilities.PairPotentials
using Unitful
using Base.Threads


hydrogen = isolated_system([
    :H => [0.1, 0, 0.]u"Å",
    :H => [0, 0, 1.]u"Å",
    :H => [4., 0, 0.]u"Å",
    :H => [4., 1., 0.]u"Å"
])

box = (
    [10.0, 0., 0.]u"Å",
    [0.0, 10., 0.]u"Å",
    [0.0, 0., 10.]u"Å",
)

pbc = (false, false, false)

hydrogen = FlexibleSystem(hydrogen[:], bounding_box=box, periodicity=pbc)


V = SimplePairPotential(
    x-> (x-0.9)^2-1,
    1,
    1,
    2.0u"Å"
)


ipi_future = @spawn IPIcalculator()
sleep(1) # we need to yeald to start the server

ipi_driver = @spawn run_driver("127.0.0.1", V, hydrogen)
sleep(1) # we need to yeald to connect to the server

calc = fetch(ipi_future)

##

test_energy_forces_virial(hydrogen, calc)

@test AtomsCalculators.potential_energy(hydrogen, V) ≈ AtomsCalculators.potential_energy(hydrogen, calc)
f_v   = AtomsCalculators.forces(hydrogen, V)
f_ipi = AtomsCalculators.forces(hydrogen, calc)
@test all( isapprox.(f_v, f_ipi) )
@test AtomsCalculators.virial(hydrogen, V) ≈ AtomsCalculators.virial(hydrogen, calc)