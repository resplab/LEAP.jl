using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_death(config)
    @testset "test Death constructor config" begin
        province = "BC"
        starting_year = 2024
        death = LEAP.Death(config["death"], province, starting_year)
        @test death.parameters[:β0] == 0
        @test death.parameters[:β1] == 0
        @test death.parameters[:β2] == 0
    end
end

function test_agent_dies(config)
    test_agent_dies_is_dead(deepcopy(config))
    test_agent_dies_is_not_dead(deepcopy(config))
end

function test_agent_dies_is_dead(config)
    @testset "test agent_dies is dead" begin
        province = "BC"
        starting_year = 2024
        config["death"]["parameters"] = Dict(
            :β0 => 1,
            :β1 => 1,
            :β2 => 1
        )
        death = LEAP.Death(config["death"], province, starting_year)
        agent = LEAP.Agent(
            sex=true,
            age=110,
            year=2024,
            year_index=1,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=1,
            SSP=config["pollution"]["SSP"]
        )
        is_dead = LEAP.agent_dies(agent, death)
        @test is_dead == true
    end
end

function test_agent_dies_is_not_dead(config)
    @testset "test agent_dies is not dead" begin
        province = "BC"
        starting_year = 2024
        config["death"]["parameters"] = Dict(
            :β0 => 0,
            :β1 => 0,
            :β2 => 0
        )
        death = LEAP.Death(config["death"], province, starting_year)
        agent = LEAP.Agent(
            sex=true,
            age=7,
            year=2025,
            year_index=2,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=1,
            SSP=config["pollution"]["SSP"]
        )
        is_dead = LEAP.agent_dies(agent, death)
        @test is_dead == false
    end
end


@testset "test death" begin
    config = JSON.parsefile(CONFIG_PATH)
    test_death(deepcopy(config))
    test_agent_dies(deepcopy(config))
end
