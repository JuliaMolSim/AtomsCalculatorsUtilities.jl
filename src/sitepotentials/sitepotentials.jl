
module SitePotentials 

using AtomsBase, AtomsCalculators, ChunkSplitters, Folds, NeighbourLists, 
      Unitful, Bumper, StrideArrays 

using StaticArrays: SVector, SMatrix 

# these have to be overloaded by new implementations of site potentials 
import AtomsCalculators: energy_unit, length_unit

# these functions will be overloaded with new methods 
import AtomsCalculators: potential_energy,
                         forces, 
                         virial, 
                         energy_forces, 
                         energy_forces_virial, 
                         force_unit 

# the following functiona are used in the implementation: 
import AtomsCalculators: zero_energy, zero_forces, zero_virial                         

include("interface.jl")

include("neighbourlist.jl")

include("assembly.jl")

include("hessian.jl")

end 
