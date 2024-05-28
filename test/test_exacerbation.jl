using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_exacerbation()
    @testset "test Exacerbation constructor config" begin
        province = "BC"
        config = JSON.parsefile(CONFIG_PATH)
        exacerbation = LEAP.Exacerbation(config["exacerbation"], province)
        @test exacerbation.hyperparameters[:β0_μ] == 0.0
        @test exacerbation.hyperparameters[:β0_σ] == 0.0000001
        @test exacerbation.parameters[:βage] == 0.0
        @test exacerbation.parameters[:βsex] == 0.0
        @test round(exacerbation.parameters[:βcontrol_C], sigdigits=4) == -1.671
        @test round(exacerbation.parameters[:βcontrol_PC], sigdigits=4) == -0.9781
        @test round(exacerbation.parameters[:βcontrol_UC], sigdigits=4) == -0.5727
        @test exacerbation.parameters[:β0] < 50
        @test exacerbation.parameters[:β0] > -50
        @test exacerbation.parameters[:β0_calibration] == 0.0
        @test exacerbation.parameters[:βasthmaDx] == 0.0
        @test exacerbation.parameters[:βprev_exac1] == 0.0
        @test exacerbation.parameters[:βprev_exac2] == 0.0
        @test exacerbation.parameters[:min_year] == 2001
        @test exacerbation.calibration_table[(2000,0)][1, :age] == 3
    end
end

function test_compute_num_exacerbations()
    @testset "test compute_num_exacerbations" begin
        province = "BC"
        age = 4
        sex = false
        year = 2001
        config = JSON.parsefile(CONFIG_PATH)
        config["exacerbation"]["hyperparameters"]["β0_μ"] = 2.0
        exacerbation = LEAP.Exacerbation(config["exacerbation"], province)
        control_levels = Dict(
            :fully_controlled => 0.0,
            :partially_controlled => 0.0,
            :uncontrolled => 1.0,
            :as_array => [0.0, 0.0, 1.0]
        )
        num_exacerbations = LEAP.compute_num_exacerbations(
            age, sex, year, control_levels, exacerbation
        )
        @test num_exacerbations > 1
    end
end

@testset "test exacerbation" begin
    test_exacerbation()
    test_compute_num_exacerbations()
end
