using Test
using LEAP
using DataFrames
using JSON

function test_incidence()
    @testset "test Incidence constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["incidence"]["hyperparameters"] = Dict(
            "β0_μ" => 0.0,
            "β0_σ" => 0.00000001
        )
        config["incidence"]["parameters"] = Dict(
            "β0" => 34.6,
            "βsex" => -9.5,
            "βage" => [-6.6, 7.7, -5.6, 3.9, -1.3],
            "βyear" => -0.019,
            "βsexage" => [-4.4, 4.7, -2.6, 0.79, 0.95],
            "βsexyear" => 0.0046,
            "βfam_hist" => [0.12, 0.36],
            "βabx_exp" => [1.8, -0.29, 0.053]
        )
        config["incidence"]["max_age"] = 60
        incidence = LEAP.Incidence(config["incidence"])
        @test incidence.hyperparameters[:β0_μ] == 0.0
        @test incidence.hyperparameters[:β0_σ] == 0.00000001
        @test incidence.parameters[:β0] == 34.6
        @test incidence.parameters[:βage] == [-6.6, 7.7, -5.6, 3.9, -1.3]
        @test incidence.parameters[:βsex] == -9.5
        @test incidence.parameters[:βyear] == -0.019
        @test incidence.parameters[:βsexage] == [-4.4, 4.7, -2.6, 0.79, 0.95]
        @test incidence.parameters[:βsexyear] == 0.0046
        @test incidence.parameters[:βfam_hist] == [0.12, 0.36]
        @test incidence.parameters[:βabx_exp] == [1.8, -0.29, 0.053]
        @test incidence.max_age == 60
    end
end


"""
    test_agent_has_asthma()

Setting the incidence parameter `βfam_hist` to [100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 1. The minimum age for
an asthma diagnosis is 3, so setting the agent age to 24 ensures that the agent is diagnosed with
asthma.

"""
function test_agent_has_asthma()
    @testset "test agent_has_asthma" begin
        max_age = 111
        province = "CA"
        starting_year = 2024
        year = 2025
        year_index = year - starting_year + 1
        month = 1
        n = 10
        population_growth_type = "M3"
        config = JSON.parsefile(LEAP.CONFIG_PATH)
        config["incidence"]["parameters"]["βfam_hist"] = [100, 0]
        config["family_history"]["parameters"]["p"] = 1.0
        incidence = LEAP.Incidence(config["incidence"])
        prevalence = LEAP.Prevalence(config["prevalence"])
        birth = LEAP.Birth(starting_year, province, population_growth_type)

        agent = LEAP.Agent(
            sex=false,
            age=24,
            year=year,
            year_index=year_index,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=month,
            SSP=config["pollution"]["SSP"]
        )
        @test agent.has_asthma == false
        @test agent_has_asthma(agent, incidence, prevalence) == true
    end
end


@testset "test incidence" begin
    test_incidence()
    test_agent_has_asthma()
end
