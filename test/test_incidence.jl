using Test
using LEAP
using DataFrames
using JSON

@testset "test agent_has_asthma" begin
    max_age = 111
    province = "CA"
    starting_year = 2024
    cal_year = 2025
    cal_year_index = cal_year - starting_year + 1
    month = 1
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.Incidence(config["incidence"])
    birth = LEAP.Birth(starting_year, province, population_growth_type)

    agent = LEAP.Agent(
        sex=false,
        age=0,
        cal_year=cal_year,
        cal_year_index=cal_year_index,
        family_hist=LEAP.FamilyHistory(config["family_history"]),
        antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
        province=province,
        month=month,
        SSP=config["pollution"]["SSP"]
    )
    @test agent.has_asthma == false
    agent_has_asthma(agent, incidence, "prevalence")
end


@testset "test agent_has_asthma error" begin
    max_age = 111
    province = "CA"
    starting_year = 2024
    cal_year = 2025
    cal_year_index = cal_year - starting_year + 1
    month = 1
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.Incidence(config["incidence"])
    birth = LEAP.Birth(starting_year, province, population_growth_type)

    agent = LEAP.Agent(
        sex=false,
        age=0,
        cal_year=cal_year,
        cal_year_index=cal_year_index,
        family_hist=LEAP.FamilyHistory(config["family_history"]),
        antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
        province=province,
        month=month,
        SSP=config["pollution"]["SSP"]
    )

    @test_throws ArgumentError(
        "inc_or_prev must be either 'prevalence' or 'incidence', received ."
    ) agent_has_asthma(agent, incidence, "")
end
