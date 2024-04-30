using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_agent()
    @testset "test Agent constructor 1 family_hist and antibiotic_exposure" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["antibiotic_exposure"]["parameters"] = Dict(
            "β0" => -100000,
            "βcal_year" => -0.01,
            "βsex" => -1,
            "θ" => 500,
            "fixyear" => nothing,
            "βfloor" => 0.0,
            "β2005" => 1,
            "β2005_cal_year" => 1
        )
        config["family_history"]["parameters"]["p"] = 0.0
        agent = LEAP.Agent(
            sex=false,
            age=23,
            cal_year=2024,
            cal_year_index=1,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"])
        )
        @test agent.sex == false
        @test agent.age == 23
        @test agent.cal_year == 2024
        @test agent.cal_year_index == 1
        @test round(agent.num_antibiotic_use, digits=1) == 0
        @test agent.has_family_hist == false
    end
    @testset "test Agent constructor 1 no family_hist or antibiotic_exposure" begin
        config = JSON.parsefile(CONFIG_PATH)
        agent = LEAP.Agent(
            sex=false,
            age=23,
            cal_year=2024,
            cal_year_index=1
        )
        @test agent.sex == false
        @test agent.age == 23
        @test agent.cal_year == 2024
        @test agent.cal_year_index == 1
        @test agent.num_antibiotic_use == nothing
        @test agent.has_family_hist == nothing
    end
end


@testset "test agent" begin
    test_agent()
end
