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


function test_compute_prob_emigration()
    @testset "test compute_prob_emigration age 0" begin
        province = "BC"
        age = 0
        sex = false
        starting_year = 2020
        year_index = 4
        population_growth_type = "FA"
        emigration = LEAP.Emigration(starting_year, province, population_growth_type)
        @test LEAP.compute_prob_emigration(year_index, age, sex, emigration) == 0
    end
    @testset "test compute_prob_emigration non-zero" begin
        province = "BC"
        age = 99
        sex = true
        starting_year = 2020
        year_index = 4
        population_growth_type = "M2"
        emigration = LEAP.Emigration(starting_year, province, population_growth_type)
        true_count = 0
        for i in 1:1000000
            if LEAP.compute_prob_emigration(year_index, age, sex, emigration)
                true_count += 1
            end
        end
        @test round(true_count, sigdigits=3) <= 436000
        @test round(true_count, sigdigits=3) >= 433000
    end
end

@testset "test emigration" begin
    test_emigration()
    test_compute_prob_emigration()
end
