using Test
using LEAP
using DataFrames
using JSON

@testset "test agent_has_asthma" begin
    max_age = 111
    province = "CA"
    starting_year = 2001
    month = 1
    time_horizon = 40
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.Incidence(config["incidence"], starting_year, province)
    birth = LEAP.Birth(starting_year, province, population_growth_type)
    antibiotic_exposure = LEAP.AntibioticExposure(config["antibiotic_exposure"])
    family_history = LEAP.FamilyHistory(config["family_history"])
    census_table = LEAP.CensusTable(config["census_table"])
    pollution_table = LEAP.PollutionTable()
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = create_agent(
        cal_year=cal_year,
        cal_year_index=tmp_cal_year_index,
        month=month,
        age=0,
        antibiotic_exposure=antibiotic_exposure,
        family_hist=family_history,
        sex=false,
        province=province,
        census_table=census_table,
        pollution_table=pollution_table
    )
    @test agent.has_asthma == false
    agent_has_asthma(agent, incidence, "prevalence")
end


@testset "test agent_has_asthma error" begin
    max_age = 111
    province = "CA"
    starting_year = 2001
    month = 1
    time_horizon = 40
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.Incidence(config["incidence"], starting_year, province)
    birth = LEAP.Birth(starting_year, province, population_growth_type)
    antibiotic_exposure = LEAP.AntibioticExposure(config["antibiotic_exposure"])
    family_history = LEAP.FamilyHistory(config["family_history"])
    census_table = LEAP.CensusTable(config["census_table"])
    pollution_table = LEAP.PollutionTable()
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = create_agent(
        cal_year=cal_year,
        cal_year_index=tmp_cal_year_index,
        month=month,
        age=0,
        antibiotic_exposure=antibiotic_exposure,
        family_hist=family_history,
        sex=false,
        province=province,
        census_table=census_table,
        pollution_table=pollution_table
    )

    @test_throws ArgumentError(
        "inc_or_prev must be either 'prevalence' or 'incidence', received ."
    ) agent_has_asthma(agent, incidence, "")
end
