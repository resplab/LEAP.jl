using Test
using LEAP
using DataFrames, JSON
using Setfield


@testset "test AsthmaCost constructor config" begin
    config = JSON.parsefile(CONFIG_PATH)
    cost = LEAP.AsthmaCost(config["cost"])
    @test cost.parameters[:control] == Array{Float64, 1}([2372, 2965, 3127]) * 1.66
    @test cost.parameters[:exac] == Array{Float64, 1}([130, 594, 2425, 9900]) * 1.66
end


@testset "test compute_cost has_asthma" begin
    config = JSON.parsefile(CONFIG_PATH)
    asthma_cost = LEAP.AsthmaCost(config["cost"])
    agent = LEAP.Agent(
        sex=true,
        age=20,
        cal_year=2021,
        cal_year_index=1,
        family_hist=LEAP.FamilyHistory(config["family_history"]),
        antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
        province="BC",
        month=1,
        SSP=config["pollution"]["SSP"]
    )
    @set agent.has_asthma = true
    cost = LEAP.compute_cost(agent, asthma_cost)
    @test cost > 0.0
end


@testset "test compute_cost !has_asthma" begin
    config = JSON.parsefile(CONFIG_PATH)
    asthma_cost = LEAP.AsthmaCost(config["cost"])
    agent = LEAP.Agent(
        sex=true,
        age=20,
        cal_year=2021,
        cal_year_index=1,
        family_hist=LEAP.FamilyHistory(config["family_history"]),
        antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
        province="BC",
        month=1,
        SSP=config["pollution"]["SSP"]
    )
    @set agent.has_asthma = false
    cost = LEAP.compute_cost(agent, asthma_cost)
    @test == 0.0
end
