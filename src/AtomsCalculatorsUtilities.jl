module AtomsCalculatorsUtilities



module Calculators
    using AtomsCalculators
    using AtomsBase
    using Unitful

    include("calculators/combination_calculator.jl")
    include("calculators/reporting_calculator.jl")
    include("calculators/subsystem_calculator.jl")
    include("calculators/zero_virial_calculator.jl")
end # module Calculators

module IPI
    include("ipi/ipi_interface.jl")
end

module Testing 
    include("testing/fdtests.jl")
end 

end
