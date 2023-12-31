using Test
using LEAP


@testset "All Tests" begin

include("test_simulation.jl")
include("test_incidence.jl")
include("test_birth.jl")

end
