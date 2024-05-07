using Test
using LEAP
using DataFrames
using Setfield


function test_birth()
    @testset "test Birth constructor no data" begin
        starting_year = 2022
        province = "BC"
        population_growth_type = "LG"
        birth = LEAP.Birth(starting_year, province, population_growth_type)
        @test birth.estimate[1, :N_relative] == 1.0
        @test birth.estimate[1, :province] == province
        @test birth.estimate[1, :N] == 42200
        @test round(birth.estimate[1, :prop_male], sigdigits=5) == 0.51185
        @test birth.estimate[1, :projection_scenario] == population_growth_type
        @test birth.estimate[1, :year] == starting_year
    end
end

function test_get_initial_population_indices()
    @testset "test get_initial_population_indices" begin
        max_age = 2
        province = "BC"
        starting_year = 2024
        population_growth_type = "M3"
        birth = LEAP.Birth(starting_year, province, population_growth_type)
        @set! birth.initial_population = DataFrames.DataFrame(
            [[0, 1, 2], [1.0, 2.0, 0.5]],
            [:age, :prop]
        )
        initial_pop_indices = get_initial_population_indices(birth, 2)
        @test initial_pop_indices == [1, 1, 2, 2, 2, 2, 3]
    end
end

function test_get_num_newborn()
    @testset "test get_num_newborn" begin
        max_age = 2
        province = "BC"
        starting_year = 2022
        cal_year_index = 3
        num_births_initial = 1000
        population_growth_type = "LG"
        birth = LEAP.Birth(starting_year, province, population_growth_type)
        num_new_born = LEAP.get_num_newborn(birth, num_births_initial, cal_year_index)
        @test num_new_born == 982
    end
end


@testset "test birth" begin
    test_birth()
    test_get_initial_population_indices()
    test_get_num_newborn()
end
