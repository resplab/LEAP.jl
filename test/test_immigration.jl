using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_immigration()
    @testset "test Immigration constructor no data" begin
        province = "BC"
        starting_year = 2024
        population_growth_type = "LG"
        max_age = 111
        immigration = LEAP.Immigration(starting_year, province, population_growth_type, max_age)
        @test round(immigration.table[(2026,)][7, :prop_immigrants_birth], sigdigits=4) == 0.004980
        @test immigration.table[(2026,)][5, :age] == ceil(Int, 5/2)
        @test round(immigration.table[(2025,)][8, :prop_immigrants_year], sigdigits=4) == 0.007294
    end
end

function test_get_num_new_immigrants()
    @testset "test get_num_new_immigrants" begin
        province = "BC"
        starting_year = 2024
        population_growth_type = "LG"
        max_age = 111
        immigration = LEAP.Immigration(starting_year, province, population_growth_type, max_age)
        year_index = 2
        num_new_born = 1000
        @test LEAP.get_num_new_immigrants(immigration, num_new_born, year_index) == 1013
    end
end

@testset "test immigration" begin
    test_immigration()
    test_get_num_new_immigrants()
end
