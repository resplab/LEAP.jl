using Test
using LEAP
using DataFrames, JSON
using Setfield


@testset "test AntibioticExposure constructor config" begin
    config = JSON.parsefile(CONFIG_PATH)
    antibiotic_exposure = LEAP.AntibioticExposure(config["antibiotic_exposure"])
    @test antibiotic_exposure.parameters[:β0] == 110.000442
    @test antibiotic_exposure.mid_trends[(2002, 0)].year[1] == 2002
end


@testset "test antibiotic_exposure_prob" begin
    parameters = Dict(
        :β0 => 1,
        :βcal_year => 0.01,
        :βsex => 1,
        :θ => 500,
        :fixyear => nothing,
        :βfloor => 0.05,
        :β2005 => 1,
        :β2005_cal_year => 1
    )
    p = LEAP.antibiotic_exposure_prob(sex=false, cal_year=2001, parameters=parameters)
    @test round(p, sigdigits=3) == 0.000000375
end


@testset "test compute_num_antibiotic_use" begin
    config = JSON.parsefile(CONFIG_PATH)
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
    antibiotic_exposure = LEAP.AntibioticExposure(config["antibiotic_exposure"])
    num_antibiotic_use = LEAP.compute_num_antibiotic_use(
        antibiotic_exposure=antibiotic_exposure,
        sex=false,
        birth_year=2001
    )
    @test round(num_antibiotic_use, digits=1) == 0
end
