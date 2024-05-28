using Test
using LEAP
using DataFrames
using JSON

function test_prevalence()
    @testset "test Prevalence constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["prevalence"]["hyperparameters"] = Dict(
            "β0_μ" => 0.0,
            "β0_σ" => 0.00000001
        )
        config["prevalence"]["parameters"] = Dict(
            "β0" => -2.28,
            "βsex" => -0.11,
            "βage" => [1.7, -2.1, 3.6, -2.9, 1.4],
            "βyear" => [2.8, -1.1],
            "βsexage" => [-7.69, 2.68, 0.86, -0.656, -0.027],
            "βsexyear" => [1.29, 0.036],
            "βyearage" => [50.6, 6.5, -39.4, 3.6, 15.9, -4.7, -7.1, 4.1, -4.8, -3.3],
            "βsexyearage" => [-3.1, 7.2, -25.7, 0.2, 11.3, -2.5, 7.6, 4.1, -15.2, 3.7],
            "βfam_hist" => [0.122, 0.376],
            "βabx_exp" => [1.826, -0.225, 0.053]
        )
        config["prevalence"]["max_age"] = 60
        prevalence = LEAP.Prevalence(config["prevalence"])
        @test prevalence.hyperparameters[:β0_μ] == 0.0
        @test prevalence.hyperparameters[:β0_σ] == 0.00000001
        @test prevalence.parameters[:β0] == -2.28
        @test prevalence.parameters[:βage] == [1.7, -2.1, 3.6, -2.9, 1.4]
        @test prevalence.parameters[:βsex] == -0.11
        @test prevalence.parameters[:βyear] == [2.8, -1.1]
        @test prevalence.parameters[:βsexage] == [-7.69, 2.68, 0.86, -0.656, -0.027]
        @test prevalence.parameters[:βsexyear] == [1.29, 0.036]
        @test prevalence.parameters[:βyearage] == [
            50.6, 6.5, -39.4, 3.6, 15.9, -4.7, -7.1, 4.1, -4.8, -3.3
        ]
        @test prevalence.parameters[:βsexyearage] == [
            -3.1, 7.2, -25.7, 0.2, 11.3, -2.5, 7.6, 4.1, -15.2, 3.7
        ]
        @test prevalence.parameters[:βfam_hist] == [0.122, 0.376]
        @test prevalence.parameters[:βabx_exp] == [1.826, -0.225, 0.053]
        @test prevalence.max_age == 60
    end
end

function test_agent_has_asthma()
    @testset "test agent_has_asthma no asthma age < 3" begin
        max_age = 111
        province = "CA"
        starting_year = 2024
        year = 2025
        year_index = year - starting_year + 1
        month = 1
        n = 10
        population_growth_type = "M3"
        config = JSON.parsefile(LEAP.CONFIG_PATH)
        prevalence = LEAP.Prevalence(config["prevalence"])
        birth = LEAP.Birth(starting_year, province, population_growth_type)

        agent = LEAP.Agent(
            sex=false,
            age=0,
            year=year,
            year_index=year_index,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=month,
            SSP=config["pollution"]["SSP"]
        )
        @test agent.has_asthma == false
        @test agent_has_asthma(agent, prevalence) == false
    end
    @testset "test agent_has_asthma asthma age >= 3" begin
        max_age = 111
        province = "CA"
        starting_year = 2024
        year = 2025
        year_index = year - starting_year + 1
        month = 1
        n = 10
        population_growth_type = "M3"
        config = JSON.parsefile(LEAP.CONFIG_PATH)
        config["prevalence"]["parameters"]["βfam_hist"] = [100, 0]
        config["family_history"]["parameters"]["p"] = 1.0
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
        @test agent_has_asthma(agent, prevalence) == true
    end
end


@testset "test prevalence" begin
    test_prevalence()
    test_agent_has_asthma()
end
