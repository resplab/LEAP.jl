using Test
using LEAP
using DataFrames, JSON, Setfield


function test_run_simulation(config)
    test_run_simulation_full(deepcopy(config))
    test_run_simulation_one_year(deepcopy(config))
end


function test_run_simulation_full(config)
    @testset "test run_simulation" begin
        config["simulation"] = Dict(
            "min_cal_year" => 2024,
            "time_horizon" => 3,
            "province" => "CA",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => 111
        )
        outcome_matrix = LEAP.run_simulation(seed=1, until_all_die=false, verbose=false, config=config)
    end
end


"""
    test_run_simulation_one_year(config)

Setting the incidence parameter `βfam_hist` to [100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 1. The maximum age is set
to 4, and the minimum age required for an asthma diagnosis is 3. So all agents aged 4 should
receive an asthma diagnosis.

Setting the exacerbation hyperparameter `β0_μ` to 20.0 ensures that every agent aged 4 has an
asthma exacerbation.

Setting the `time_horizon` to 1 means that the agents are generated from the initial population
table, and that no immigration happens.

Setting the antibiotic exposure parameters below ensures that the antibiotic use is 0.

Setting the `num_births_initial` to 10 and starting in 2024 with growth type "M3", each of the
age groups has 10 agents, for a total of 10 x 5 = 50 agents.

Setting the `prevalence` parameters below ensures that the prevalence is 0.

"""
function test_run_simulation_one_year(config)
    @testset "test run_simulation one year" begin
        config["simulation"] = Dict(
            "min_cal_year" => 2024,
            "time_horizon" => 1,
            "province" => "CA",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => 4
        )
        config["antibiotic_exposure"]["parameters"] = Dict(
            :β0 => -100000,
            :βcal_year => -0.01,
            :βsex => -1,
            :θ => 500,
            :fixyear => nothing,
            :βfloor => 0.0,
            :β2005 => 1,
            :β2005_cal_year => 1
        )
        config["exacerbation"]["hyperparameters"]["β0_μ"] = 5.0
        config["incidence"]["parameters"]["βfam_hist"] = [100, 0]
        config["prevalence"]["parameters"] = Dict(
            "β0" => -20,
            "βsex" => -20,
            "βage" => [0.0, 0.0, 0.0, 0.0, 0.0],
            "βyear" => [0.0, 0.0],
            "βsexage" => [0.0, 0.0, 0.0, 0.0, 0.0],
            "βsexyear" => [0.0, 0.0],
            "βyearage" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "βsexyearage" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "βfam_hist" => [-100, 0],
            "βabx_exp" => [0.0, 0.0, 0.0]
        )
        config["family_history"]["parameters"]["p"] = 1.0
        config["death"]["parameters"] = Dict(
            :β0 => -1,
            :β1 => -1,
            :β2 => -1
        )
        outcome_matrix = LEAP.run_simulation(seed=1, until_all_die=false, verbose=true, config=config)
        @info "asthma_incidence: $(outcome_matrix.asthma_incidence)"
        @test outcome_matrix.immigration == zeros(Int, (1, 5, 2))
        @test outcome_matrix.antibiotic_exposure == zeros(Int, (1, 5, 2))
        @test outcome_matrix.death == zeros(Int, (1, 5, 2))
        @test (
            outcome_matrix.asthma_incidence[1, 5, 1] +
            outcome_matrix.asthma_incidence[1, 5, 2]
        ) == 10
        @test (
            outcome_matrix.exacerbation[1, 5, 1] +
            outcome_matrix.exacerbation[1, 5, 2]
        ) > 1000
    end
end


"""
    test_check_if_agent_gets_new_asthma_diagnosis!(config)

Setting the `time_horizon` to 1 means that the agents are generated from the initial population
table, and that no immigration happens.

Setting the antibiotic exposure parameters below ensures that the antibiotic use is 0.

Setting the `num_births_initial` to 10 and starting in 2024 with growth type "M3", each of the
age groups has 10 agents, for a total of 10 x 5 = 50 agents.

testset 1:
==========
Setting the incidence parameter `βfam_hist` to [100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 1. The maximum age is set
to 4, and the minimum age required for an asthma diagnosis is 3. So all agents aged 4 should
receive an asthma diagnosis.

testset 2:
==========
Setting the incidence parameter `βfam_hist` to [-100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 0.

"""
function test_check_if_agent_gets_new_asthma_diagnosis!(config)
    @testset "test check_if_agent_gets_new_asthma_diagnosis! true" begin
        min_cal_year = 2024
        time_horizon = 1
        max_cal_year = min_cal_year + time_horizon - 1
        max_age = 4

        config["simulation"] = Dict(
            "min_cal_year" => min_cal_year,
            "time_horizon" => time_horizon,
            "province" => "CA",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => max_age
        )
        config["antibiotic_exposure"]["parameters"] = Dict(
            :β0 => -100000,
            :βcal_year => -0.01,
            :βsex => -1,
            :θ => 500,
            :fixyear => nothing,
            :βfloor => 0.0,
            :β2005 => 1,
            :β2005_cal_year => 1
        )
        config["incidence"]["parameters"]["βfam_hist"] = [100, 0]
        config["family_history"]["parameters"]["p"] = 1.0
        outcome_matrix = LEAP.create_outcome_matrix(
            until_all_die=false,
            cal_years=min_cal_year:max_cal_year,
            min_cal_year=min_cal_year,
            max_cal_year=max_cal_year,
            max_age=max_age
        )
        simulation = LEAP.Simulation(config)
        @set! simulation.agent = LEAP.Agent(
            sex=false,
            age=4,
            cal_year=min_cal_year,
            cal_year_index=1,
            family_hist=simulation.family_history,
            antibiotic_exposure=simulation.antibiotic_exposure,
            province=simulation.province,
            SSP=simulation.SSP
        )
        @test simulation.agent !== nothing
        LEAP.check_if_agent_gets_new_asthma_diagnosis!(simulation, outcome_matrix)
        @info "asthma_incidence: $(outcome_matrix.asthma_incidence)"
        @test outcome_matrix.asthma_incidence[1, 5, 1] == 1
        @test outcome_matrix.asthma_incidence[1, 5, 2] == 0
        @test outcome_matrix.asthma_status[1, 5, 1] == 1
        @test outcome_matrix.asthma_status[1, 5, 2] == 0
        @test simulation.agent.has_asthma == true
        @test simulation.agent.asthma_age == 4
        @test simulation.agent.asthma_status == true
    end
    @testset "test check_if_agent_gets_new_asthma_diagnosis! false" begin
        min_cal_year = 2024
        time_horizon = 1
        max_cal_year = min_cal_year + time_horizon - 1
        max_age = 4

        config["simulation"] = Dict(
            "min_cal_year" => min_cal_year,
            "time_horizon" => time_horizon,
            "province" => "CA",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => max_age
        )
        config["antibiotic_exposure"]["parameters"] = Dict(
            :β0 => -100000,
            :βcal_year => -0.01,
            :βsex => -1,
            :θ => 500,
            :fixyear => nothing,
            :βfloor => 0.0,
            :β2005 => 1,
            :β2005_cal_year => 1
        )
        config["incidence"]["parameters"] = Dict(
            "β0" => 34.6,
            "βsex" => -9.5,
            "βage" => [-6.6, 7.7, -5.6, 3.9, -1.3],
            "βyear" => -0.019,
            "βsexage" => [-4.4, 4.7, -2.6, 0.79, 0.95],
            "βsexyear" => 0.0046,
            "βfam_hist" => [-100, 0],
            "βabx_exp" => [0.0, 0.0, 0.0]
        )
        config["family_history"]["parameters"]["p"] = 1.0
        outcome_matrix = LEAP.create_outcome_matrix(
            until_all_die=false,
            cal_years=min_cal_year:max_cal_year,
            min_cal_year=min_cal_year,
            max_cal_year=max_cal_year,
            max_age=max_age
        )
        simulation = LEAP.Simulation(config)
        @set! simulation.agent = LEAP.Agent(
            sex=false,
            age=4,
            cal_year=min_cal_year,
            cal_year_index=1,
            family_hist=simulation.family_history,
            antibiotic_exposure=simulation.antibiotic_exposure,
            province=simulation.province,
            SSP=simulation.SSP
        )
        @test simulation.agent !== nothing
        LEAP.check_if_agent_gets_new_asthma_diagnosis!(simulation, outcome_matrix)
        @info "asthma_incidence: $(outcome_matrix.asthma_incidence)"
        @test outcome_matrix.asthma_incidence[1, 5, 1] == 0
        @test outcome_matrix.asthma_incidence[1, 5, 2] == 0
        @test outcome_matrix.asthma_status[1, 5, 1] == 0
        @test outcome_matrix.asthma_status[1, 5, 2] == 0
        @test simulation.agent.has_asthma == false
        @test simulation.agent.asthma_age == nothing
        @test simulation.agent.asthma_status == false
    end
end


"""
    test_update_asthma_effects!(config)

Setting the `time_horizon` to 1 means that the agents are generated from the initial population
table, and that no immigration happens.

Setting the antibiotic exposure parameters below ensures that the antibiotic use is 0.

Setting the `num_births_initial` to 10 and starting in 2024 with growth type "M3", each of the
age groups has 10 agents, for a total of 10 x 5 = 50 agents.

Setting the incidence parameter `βfam_hist` to [100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 1. The maximum age is set
to 4, and the minimum age required for an asthma diagnosis is 3. So all agents aged 4 should
receive an asthma diagnosis.

Setting the control parameter `θ` to [-1e5, -1e5] ensures that the `control_levels` are:
    FC: 0.0
    PC: 0.0
    UC: 1.0

Setting the exacerbation parameter `β0_μ` to 5.0 ensures that the number of exacerbations will be
large.

"""
function test_update_asthma_effects!(config)
    @testset "test update_asthma_effects!" begin
        min_cal_year = 2024
        time_horizon = 1
        max_cal_year = min_cal_year + time_horizon - 1
        max_age = 4

        config["simulation"] = Dict(
            "min_cal_year" => min_cal_year,
            "time_horizon" => time_horizon,
            "province" => "CA",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => max_age
        )
        config["antibiotic_exposure"]["parameters"] = Dict(
            :β0 => -100000,
            :βcal_year => -0.01,
            :βsex => -1,
            :θ => 500,
            :fixyear => nothing,
            :βfloor => 0.0,
            :β2005 => 1,
            :β2005_cal_year => 1
        )
        config["incidence"]["parameters"]["βfam_hist"] = [100, 0]
        config["family_history"]["parameters"]["p"] = 1.0
        config["control"]["parameters"]["θ"] = [-1e5, -1e5]
        config["exacerbation"]["hyperparameters"]["β0_μ"] = 5.0
        outcome_matrix = LEAP.create_outcome_matrix(
            until_all_die=false,
            cal_years=min_cal_year:max_cal_year,
            min_cal_year=min_cal_year,
            max_cal_year=max_cal_year,
            max_age=max_age
        )
        simulation = LEAP.Simulation(config)
        @set! simulation.agent = LEAP.Agent(
            sex=false,
            age=4,
            cal_year=min_cal_year,
            cal_year_index=1,
            family_hist=simulation.family_history,
            antibiotic_exposure=simulation.antibiotic_exposure,
            province=simulation.province,
            SSP=simulation.SSP,
            has_asthma=true,
            asthma_age=4,
            asthma_status=true
        )
        @test simulation.agent !== nothing
        LEAP.update_asthma_effects!(simulation, outcome_matrix)
        @info "control levels: $(outcome_matrix.control)"
        @test simulation.agent.has_asthma == true
        @test simulation.agent.asthma_age == 4
        @test simulation.agent.asthma_status == true
        @test simulation.agent.control_levels == Dict(
            :fully_controlled => 0.0,
            :partially_controlled => 0.0,
            :uncontrolled => 1.0,
            :as_array => [0.0, 0.0, 1.0]
        )
        @test outcome_matrix.control[1, 5, 1, :] == [0.0, 0.0, 1.0]
        @test simulation.agent.exac_hist.num_current_year > 100
        @test simulation.agent.exac_hist.num_prev_year == 0
        @test outcome_matrix.exacerbation[1, 5, 1] > 100
    end
end


"""
    test_reassess_asthma_diagnosis!(config)

Setting the `time_horizon` to 1 means that the agents are generated from the initial population
table, and that no immigration happens.

Setting the antibiotic exposure parameters below ensures that the antibiotic use is 0.

Setting the `num_births_initial` to 10 and starting in 2024 with growth type "M3", each of the
age groups has 10 agents, for a total of 10 x 5 = 50 agents.

Setting the incidence parameter `βfam_hist` to [100, 0] and the family history parameter `p` to 1.0
ensures that the probability of an agent being diagnosed with asthma is 1. The maximum age is set
to 4, and the minimum age required for an asthma diagnosis is 3. So all agents aged 4 should
receive an asthma diagnosis.

Setting the control parameter `θ` to [-1e5, -1e5] ensures that the `control_levels` are:
    FC: 0.0
    PC: 0.0
    UC: 1.0

Setting the exacerbation parameter `β0_μ` to 5.0 ensures that the number of exacerbations will be
large.

"""
function test_reassess_asthma_diagnosis!(config)
    @testset "test test_reassess_asthma_diagnosis!" begin
        min_cal_year = 2024
        time_horizon = 1
        max_cal_year = min_cal_year + time_horizon - 1
        max_age = 100

        config["simulation"] = Dict(
            "min_cal_year" => min_cal_year,
            "time_horizon" => time_horizon,
            "province" => "BC",
            "population_growth_type" => "M3",
            "num_births_initial" => 10,
            "max_age" => max_age
        )
        config["antibiotic_exposure"]["parameters"] = Dict(
            :β0 => -100000,
            :βcal_year => -0.01,
            :βsex => -1,
            :θ => 500,
            :fixyear => nothing,
            :βfloor => 0.0,
            :β2005 => 1,
            :β2005_cal_year => 1
        )
        config["incidence"]["parameters"]["βfam_hist"] = [100, 0]
        config["family_history"]["parameters"]["p"] = 1.0
        config["control"]["parameters"]["θ"] = [-1e5, -1e5]
        config["exacerbation"]["hyperparameters"]["β0_μ"] = 10.0
        outcome_matrix = LEAP.create_outcome_matrix(
            until_all_die=false,
            cal_years=min_cal_year:max_cal_year,
            min_cal_year=min_cal_year,
            max_cal_year=max_cal_year,
            max_age=max_age
        )
        simulation = LEAP.Simulation(config)
        @set! simulation.agent = LEAP.Agent(
            sex=true,
            age=53,
            cal_year=min_cal_year,
            cal_year_index=1,
            family_hist=simulation.family_history,
            antibiotic_exposure=simulation.antibiotic_exposure,
            province=simulation.province,
            SSP=simulation.SSP,
            has_asthma=true,
            asthma_age=4,
            asthma_status=true,
            exac_hist=LEAP.ExacerbationHist(20, 0)
        )
        LEAP.reassess_asthma_diagnosis!(simulation, outcome_matrix)
        @test simulation.agent.has_asthma == true
        @test simulation.agent.asthma_age == 4
        @test simulation.agent.asthma_status == true
        @test simulation.agent.control_levels == Dict(
            :fully_controlled => 0.0,
            :partially_controlled => 0.0,
            :uncontrolled => 1.0,
            :as_array => [0.0, 0.0, 1.0]
        )
        @test outcome_matrix.control[1, 54, 2, :] == [0.0, 0.0, 1.0]
        @test simulation.agent.exac_hist.num_current_year > 100
        @test simulation.agent.exac_hist.num_prev_year == 20
        @test outcome_matrix.exacerbation[1, 54, 2] > 100
    end
end


function test_get_new_agents(config)
    config["simulation"] = Dict(
        "min_cal_year" => 2024,
        "time_horizon" => 3,
        "province" => "CA",
        "population_growth_type" => "M3",
        "num_births_initial" => 10,
        "max_age" => 111
    )
    test_get_new_agents_initial_year(config)
    test_get_new_agents_subsequent_year(config)
end

function test_get_new_agents_initial_year(config)
    @testset "test get_new_agents initial cal_year" begin

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
end


function test_get_new_agents_subsequent_year(config)
    @testset "test get_new_agents subsequent cal_year" begin
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
end


@testset "test simulation" begin
    config = JSON.parsefile(CONFIG_PATH)
    # test_get_new_agents(deepcopy(config))
    # test_check_if_agent_gets_new_asthma_diagnosis!(deepcopy(config))
    test_run_simulation_one_year(deepcopy(config))
    test_update_asthma_effects!(deepcopy(config))
end
