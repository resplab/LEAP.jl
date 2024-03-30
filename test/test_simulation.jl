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
