using Test
using LEAP
using DataFrames
using Setfield

@testset "test get_initial_population_indices" begin
    max_age = 2
    province = "BC"
    starting_year = 2001
    population_growth_type = "M3"
    birth = LEAP.set_up_birth(starting_year, population_growth_type, province)
    @set! birth.initial_population = DataFrames.DataFrame(
        [[0, 1, 2], [1.0, 2.0, 0.5]],
        [:age, :prop]
    )
    initial_pop_indices = get_initial_population_indices(birth, 2)
    @test initial_pop_indices == [1, 1, 2, 2, 2, 2, 3]
end
