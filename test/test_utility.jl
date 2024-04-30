using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_utility()
    @testset "test Utility constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["utility"]["parameters"] = Dict(
            "βcontrol" => [0.06, 0.09, 0.10],
            "βexac_sev_hist" => [0.00615, 0.0169, 0.0192, 0.0215]
        )
        utility = LEAP.Utility(config["utility"])
        @test utility.parameters[:βcontrol] == [0.06, 0.09, 0.10]
        @test utility.parameters[:βexac_sev_hist] == [0.00615, 0.0169, 0.0192, 0.0215]
        @test utility.table[(9, 0)].eq5d[1] == 0.9405
    end
end

function test_compute_utility()
    @testset "test compute_utility no asthma" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["utility"]["parameters"] = Dict(
            "βcontrol" => [0.06, 0.09, 0.10],
            "βexac_sev_hist" => [0.00615, 0.0169, 0.0192, 0.0215]
        )
        utility = LEAP.Utility(config["utility"])
        agent = LEAP.Agent(
            sex=false,
            age=9,
            cal_year=2024,
            cal_year_index=1,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province="BC",
            month=1,
            SSP=config["pollution"]["SSP"],
            has_asthma=false
        )
        @test LEAP.compute_utility(agent, utility) == 0.9405
    end
    @testset "test compute_utility has asthma" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["utility"]["parameters"] = Dict(
            "βcontrol" => [0.06, 0.09, 0.0405],
            "βexac_sev_hist" => [0.0, 0.02, 0.0, 0.0]
        )
        utility = LEAP.Utility(config["utility"])
        agent = LEAP.Agent(
            sex=false,
            age=9,
            cal_year=2024,
            cal_year_index=1,
            exac_sev_hist=LEAP.ExacerbationSeverityHist([1, 5, 1, 0], zeros(4)),
            control_levels=Dict(
                :fully_controlled => 0.0,
                :partially_controlled => 0.0,
                :uncontrolled => 1.0,
                :as_array => [0.0, 0.0, 1.0]
            ),
            province="BC",
            month=1,
            SSP=config["pollution"]["SSP"],
            has_asthma=true,
            asthma_age=7
        )
        @test LEAP.compute_utility(agent, utility) == 0.8
    end
end

@testset "test utility" begin
    test_utility()
    test_compute_utility()
end
