using Test
using LEAP
using DataFrames, JSON


@testset "test simulation" begin
    config = JSON.parsefile(CONFIG_PATH)
    config["simulation"] = Dict(
        "min_cal_year" => 2024,
        "time_horizon" => 3,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    LEAP.run_simulation(seed=1, until_all_die=false, verbose=false, config=config)
end


@testset "test get_new_agents initial cal_year" begin
    config = JSON.parsefile(CONFIG_PATH)
    config["simulation"] = Dict(
        "min_cal_year" => 2024,
        "time_horizon" => 3,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    simulation = LEAP.Simulation(config)
    cal_year = 2024
    cal_year_index = cal_year - config["simulation"]["min_cal_year"] + 1
    new_agents_df = LEAP.get_new_agents(
        simulation=simulation,
        cal_year=cal_year,
        cal_year_index=cal_year_index
    )
    @test size(new_agents_df)[1] == 999
    @test size(new_agents_df)[2] == 3
    @test new_agents_df.age[1] == 0
    @test new_agents_df.age[10] == 0
    @test new_agents_df.immigrant == falses(999)
end


@testset "test get_new_agents subsequent cal_year" begin
    config = JSON.parsefile(CONFIG_PATH)
    config["simulation"] = Dict(
        "min_cal_year" => 2024,
        "time_horizon" => 3,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    simulation = LEAP.Simulation(config)
    cal_year = 2025
    cal_year_index = cal_year - config["simulation"]["min_cal_year"] + 1
    new_agents_df = LEAP.get_new_agents(
        simulation=simulation,
        cal_year=cal_year,
        cal_year_index=cal_year_index
    )
    @test size(new_agents_df)[1] == 21
    @test size(new_agents_df)[2] == 3
    @test new_agents_df.age[11] == 0
    @test new_agents_df.age[21] == 0
    @test new_agents_df.immigrant == vcat(trues(10), falses(11))
end
