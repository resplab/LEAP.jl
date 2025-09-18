using Test
using LEAP
using DataFrames
using JSON

function test_exacerbation_severity()
    @testset "test ExacerbationSeverity constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["exacerbation_severity"]["hyperparameters"] = Dict(
            "α" => [1000, 0.00001, 0.00001, 0.00001],
            "k" => 100
        )
        config["exacerbation_severity"]["parameters"] = Dict(
            "p" => [0.25, 0.25, 0.25, 0.25],
            "βprev_hosp_ped" => 1.79,
            "βprev_hosp_adult" => 2.88
        )
        exacerbation_severity = LEAP.ExacerbationSeverity(config["exacerbation_severity"])
        @test exacerbation_severity.hyperparameters[:α] == [1000, 0.00001, 0.00001, 0.00001]
        @test exacerbation_severity.hyperparameters[:k] == 100
        @test round.(exacerbation_severity.parameters[:p], digits=1) == [1.0, 0.0, 0.0, 0.0]
        @test exacerbation_severity.parameters[:βprev_hosp_ped] == 1.79
        @test exacerbation_severity.parameters[:βprev_hosp_adult] == 2.88
    end
end


"""
    test_compute_hospitalization_prob()

Setting the exacerbation hyperparameter `β0_μ` to 20.0 ensures that every agent aged 4+ has an
asthma exacerbation.

Setting the exacerbation_severty hyperparameter `α` to [0.00001, 0.00001, 0.00001, 1000] ensures
that probability of having an exacerbation with high severity is 1.

"""
function test_compute_hospitalization_prob()
    @testset "test compute_hospitalization_prob" begin
        year = 2024
        year_index = 1
        province = "BC"
        config = JSON.parsefile(CONFIG_PATH)
        config["exacerbation_severity"]["hyperparameters"] = Dict(
            "α" => [0.00001, 0.00001, 0.00001, 1000],
            "k" => 100
        )
        config["exacerbation_severity"]["parameters"] = Dict(
            "p" => [0.25, 0.25, 0.25, 0.25],
            "βprev_hosp_ped" => 1.79,
            "βprev_hosp_adult" => 2.88
        )
        config["exacerbation"]["hyperparameters"]["β0_μ"] = 20.0
        exacerbation_severity = LEAP.ExacerbationSeverity(config["exacerbation_severity"])
        control = LEAP.Control(config["control"])
        exacerbation = LEAP.Exacerbation(config["exacerbation"], province)
        agent = LEAP.Agent(
            sex=true,
            age=90,
            year=year,
            year_index=year_index,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=1,
            SSP=config["pollution"]["SSP"],
            has_asthma=true,
            asthma_age=85
        )
        prob = LEAP.compute_hospitalization_prob(
            agent, exacerbation_severity, control, exacerbation
        )
        @test prob == 1.0
    end
end



@testset "test severity" begin
    test_exacerbation_severity()
    test_compute_hospitalization_prob()
end
