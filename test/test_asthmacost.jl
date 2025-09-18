using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_asthma_cost()
    @testset "test AsthmaCost constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        cost = LEAP.AsthmaCost(config["cost"])
        @test cost.parameters[:control] == Array{Float64, 1}([2372, 2965, 3127]) * 1.66
        @test cost.parameters[:exac] == Array{Float64, 1}([130, 594, 2425, 9900]) * 1.66
    end
end

function test_compute_cost()
    @testset "test compute_cost" begin
        config = JSON.parsefile(CONFIG_PATH)
        asthma_cost = LEAP.AsthmaCost(config["cost"])
        agent = LEAP.Agent(
            sex=true,
            age=20,
            year=2024,
            year_index=1,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province="BC",
            month=1,
            SSP=config["pollution"]["SSP"]
        )
        test_compute_cost_has_asthma(deepcopy(asthma_cost), deepcopy(agent))
        test_compute_cost_no_asthma(deepcopy(asthma_cost), deepcopy(agent))
        test_compute_cost_has_asthma_missing_control_levels(asthma_cost, agent)
    end
end

function test_compute_cost_has_asthma(asthma_cost, agent)
    @testset "test compute_cost has_asthma" begin
        @set! agent.has_asthma = true
        @set! agent.control_levels = Dict(
            :fully_controlled => 0.2,
            :partially_controlled => 0.75,
            :uncontrolled => 0.05,
            :as_array => [0.2, 0.75, 0.05]
        )
        @set! agent.exac_hist.num_current_year = 1
        cost = LEAP.compute_cost(agent, asthma_cost)
        @test cost > 0.0
    end
end

function test_compute_cost_has_asthma_missing_control_levels(asthma_cost, agent)
    @testset "test compute_cost has_asthma" begin
        @set! agent.has_asthma = true
        @set! agent.exac_hist.num_current_year = 1
        control_levels = Dict(
            :fully_controlled => 0.3333,
            :partially_controlled => 0.3333,
            :uncontrolled => 0.3333,
            :as_array => [0.3333, 0.3333, 0.3333]
        )
        @test (
            @test_logs(
                (:warn, "Control levels have not been set, default to $control_levels"),
                LEAP.compute_cost(agent, asthma_cost)
        )) > 0.0
    end
end

function test_compute_cost_no_asthma(asthma_cost, agent)
    @testset "test compute_cost !has_asthma" begin
        @set! agent.has_asthma = false
        cost = LEAP.compute_cost(agent, asthma_cost)
        @test cost == 0.0
    end
end


@testset "test asthmacost" begin
    test_asthma_cost()
    test_compute_cost()
end
