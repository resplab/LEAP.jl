using Test
using LEAP
using DataFrames, JSON

@testset "test create_outcome_matrix" begin
    config = JSON.parsefile(CONFIG_PATH)
    config["simulation"] = Dict(
        "starting_year" => 2001,
        "time_horizon" => 40,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    simulation = LEAP.Simulation(config)
    until_all_die = false
    max_age = simulation.max_age
    min_cal_year = simulation.starting_calendar_year
    max_cal_year = min_cal_year + simulation.time_horizon - 1
    cal_years = min_cal_year:max_cal_year
    outcome_matrix = LEAP.create_outcome_matrix(until_all_die, cal_years, min_cal_year,
       max_cal_year, max_age)

    @test size(outcome_matrix.control) == (40, 112, 2, 3)
    @test size(outcome_matrix.exacerbation_by_severity) == (40, 112, 2, 4)
    @test typeof(outcome_matrix.asthma_prevalence_contingency_table) == GroupedDataFrame{DataFrame}
end


@testset "test add_antibiotic_use_to_outcome_matrix" begin
    config = JSON.parsefile(CONFIG_PATH)
    config["simulation"] = Dict(
        "starting_year" => 2001,
        "time_horizon" => 40,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    simulation = LEAP.Simulation(config)
    until_all_die = false
    max_age = simulation.max_age
    min_cal_year = simulation.starting_calendar_year
    max_cal_year = min_cal_year + simulation.time_horizon - 1
    cal_years = min_cal_year:max_cal_year
    outcome_matrix = LEAP.create_outcome_matrix(until_all_die, cal_years, min_cal_year,
       max_cal_year, max_age)

    LEAP.increment_field_in_outcome_matrix!(outcome_matrix, "antibiotic_exposure",
        2, false, 1, 3)
    @test outcome_matrix.antibiotic_exposure[1, 3, 1, :] == [3]
end
