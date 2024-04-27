using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_immigration()
    @testset "test Immigration constructor no data" begin
        province = "BC"
        starting_year = 2024
        population_growth_type = "LG"
        immigration = LEAP.Immigration(starting_year, province, population_growth_type)
        @test round(immigration.table[(2026,)][7, :n_prop_birth], sigdigits=4) == 0.004980
        @test immigration.table[(2026,)][5, :age] == ceil(Int, 5/2)
        @test round(immigration.table[(2025,)][8, :weights], sigdigits=4) == 0.007294
    end
end

function test_get_num_new_immigrants()
    @testset "test get_num_new_immigrants" begin
        province = "BC"
        starting_year = 2024
        population_growth_type = "LG"
        immigration = LEAP.Immigration(starting_year, province, population_growth_type)
        cal_year_index = 2
        num_new_born = 1000
        @test LEAP.get_num_new_immigrants(immigration, num_new_born, cal_year_index) == 1009
    end
end

@testset "test immigration" begin
    test_immigration()
    test_get_num_new_immigrants()
end
