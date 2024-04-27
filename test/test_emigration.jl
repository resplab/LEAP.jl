using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_emigration()
    @testset "test Emigration constructor no data" begin
        province = "BC"
        starting_year = 2024
        population_growth_type = "LG"
        emigration = LEAP.Emigration(starting_year, province, population_growth_type)
        @test emigration.table[(2025,)][1, :F] == 0.0
        @test emigration.table[(2025,)][1, :M] == 0.0
        @test round(emigration.table[(2025,)][89, :F], sigdigits=4) == 0.02836
    end
end

@testset "test emigration" begin
    test_emigration()
end
