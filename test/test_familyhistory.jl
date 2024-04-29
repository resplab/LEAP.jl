using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_familyhistory()
    @testset "test FamilyHistory constructor config" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["family_history"]["parameters"]["p"] = 0.0
        family_history = LEAP.FamilyHistory(config["family_history"])
        @test family_history.parameters[:p] == 0.0
    end
    @testset "test FamilyHistory constructor config invalid p" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["family_history"]["parameters"]["p"] = 5.6
        @test_throws ArgumentError(
            "p must be a probability between 0 and 1, received 5.6."
        ) LEAP.FamilyHistory(config["family_history"])
    end
end

function test_has_family_history_of_asthma()
    @testset "test has_family_history_of_asthma" begin
        config = JSON.parsefile(CONFIG_PATH)
        config["family_history"]["parameters"]["p"] = 0.0
        family_history = LEAP.FamilyHistory(config["family_history"])
        has_family_hist = LEAP.has_family_history_of_asthma(
            family_history
        )
        @test has_family_hist == false
    end
end

@testset "test familyhistory" begin
    test_familyhistory()
    test_has_family_history_of_asthma()
end
