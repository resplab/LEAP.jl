using Test
using LEAP


@testset "All Tests" begin

include("test_simulation.jl")
include("test_incidence.jl")
include("test_birth.jl")
include("test_censusdivision.jl")
include("test_outcomematrix.jl")
include("test_pollution.jl")

end
