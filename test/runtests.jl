using Test
using LEAP


@testset "All Tests" begin

include("test_antibioticexposure.jl")
include("test_asthmacost.jl")
include("test_birth.jl")
include("test_censusdivision.jl")
include("test_control.jl")
include("test_incidence.jl")
include("test_outcomematrix.jl")
include("test_pollution.jl")
include("test_simulation.jl")

end
