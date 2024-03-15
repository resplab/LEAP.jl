using Test
using LEAP
using DataFrames
using JSON

@testset "test agent_has_asthma" begin
    max_age = 111
    province = "CA"
    starting_year = 2001
    time_horizon = 40
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.set_up_incidence(starting_year, province)
    birth = LEAP.set_up_birth(starting_year, population_growth_type, province)
    antibiotic_exposure = LEAP.AntibioticExposure(
        config["antibiotic_exposure"], LEAP.abx_mid_trends, nothing
    )
    family_history = LEAP.set_up_family_history()
    census_table = LEAP.set_up_census_table()
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = create_agent(
        cal_year=cal_year,
        cal_year_index=tmp_cal_year_index,
        age=0,
        antibiotic_exposure=antibiotic_exposure,
        family_hist=family_history,
        sex=false,
        province=province,
        census_table=census_table
    )
    @test agent.has_asthma == false
    agent_has_asthma(agent, incidence, "prevalence")
end


@testset "test agent_has_asthma error" begin
    max_age = 111
    province = "CA"
    starting_year = 2001
    time_horizon = 40
    n = 10
    population_growth_type = "M3"
    config = JSON.parsefile(LEAP.CONFIG_PATH)
    incidence = LEAP.set_up_incidence(starting_year, province)
    birth = LEAP.set_up_birth(starting_year, population_growth_type, province)
    antibiotic_exposure = LEAP.AntibioticExposure(
        config["antibiotic_exposure"], LEAP.abx_mid_trends, nothing
    )
    family_history = LEAP.set_up_family_history()
    census_table = LEAP.set_up_census_table()
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = create_agent(
        cal_year=cal_year,
        cal_year_index=tmp_cal_year_index,
        age=0,
        antibiotic_exposure=antibiotic_exposure,
        family_hist=family_history,
        sex=false,
        province=province,
        census_table=census_table
    )

    @test_throws ArgumentError(
        "inc_or_prev must be either 'prevalence' or 'incidence', received ."
    ) agent_has_asthma(agent, incidence, "")
end
