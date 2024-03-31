using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_control()
    @testset "test Control constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        control = LEAP.Control(config["control"])
        @test control.hyperparameters[:β0_μ] == 0.0
        @test control.hyperparameters[:β0_σ] == 1.678728
        @test control.parameters[:βage] == 3.5430381
        @test control.parameters[:βage2] == -3.4980710
        @test control.parameters[:βsexage] == -0.8161495
        @test control.parameters[:βsexage2] == -1.1654264
        @test control.parameters[:βsex] == 0.2347807
        @test control.parameters[:θ] == [-0.3950, 2.754]
        @test control.parameters[:βDx2] == nothing
        @test control.parameters[:βDx3] == nothing
        @test control.parameters[:β0] > -50
        @test control.parameters[:β0] < 50
    end
end

function test_compute_control_levels()
    config = JSON.parsefile(CONFIG_PATH)
    control = LEAP.Control(config["control"])
    sex = false
    age = 20
    control_levels = LEAP.compute_control_levels(control, sex, age, false)
    @test round(sum(
        control_levels[:fully_controlled],
        control_levels[:partially_controlled],
        control_levels[:uncontrolled]
        ), digits=0) == 1
    @test round(sum(control_levels[:as_array]), digits=0) == 1
end

@testset "test control" begin
    test_control()
end
