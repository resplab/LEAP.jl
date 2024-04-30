using Test
using LEAP


@testset "All Tests" begin

include("test_antibioticexposure.jl")
include("test_asthmacost.jl")
include("test_birth.jl")
include("test_censusdivision.jl")
include("test_control.jl")
include("test_death.jl")
include("test_emigration.jl")
include("test_exacerbation.jl")
include("test_familyhistory.jl")
include("test_immigration.jl")
include("test_incidence.jl")
include("test_outcomematrix.jl")
include("test_pollution.jl")
include("test_prevalence.jl")
include("test_reassessment.jl")
include("test_severity.jl")
include("test_simulation.jl")
include("test_utility.jl")

end
