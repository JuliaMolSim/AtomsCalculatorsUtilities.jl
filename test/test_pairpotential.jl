using AtomsCalculatorsUtilities.PairPotentials
using AtomsCalculatorsUtilities.SitePotentials
using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using StaticArrays
using Test
using Unitful



hydrogen = isolated_system([
    :H => [0, 0, 0.]u"Å",
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


test_energy_forces_virial(hydrogen, V)