using Test
using LEAP
using DataFrames, JSON
using Setfield


function test_reassessment()
    @testset "test Reassessment constructor no data" begin
        starting_year = 2020
        province = "BC"
        reassessment = LEAP.Reassessment(starting_year, province)
        @test reassessment.table[1][1, :age] == 4
    end
end


"""
    test_agent_has_asthma()

By the `master_asthma_reassessment.csv` table, the probability of still having asthma is 1 for
a male aged 53 in 2024 in BC.

"""
function test_agent_has_asthma()
    @testset "test agent_has_asthma" begin
        province = "BC"
        starting_year = 2024
        cal_year = 2024
        cal_year_index = cal_year - starting_year + 1
        month = 1
        config = JSON.parsefile(LEAP.CONFIG_PATH)
        reassessment = LEAP.Reassessment(starting_year, province)

        agent = LEAP.Agent(
            sex=true,
            age=53,
            cal_year=cal_year,
            cal_year_index=cal_year_index,
            family_hist=LEAP.FamilyHistory(config["family_history"]),
            antibiotic_exposure=LEAP.AntibioticExposure(config["antibiotic_exposure"]),
            province=province,
            month=month,
            SSP=config["pollution"]["SSP"],
            has_asthma=true
        )
        @test agent.has_asthma == true
        @test agent_has_asthma(agent, reassessment) == true
    end
end

@testset "test reassessment" begin
    test_reassessment()
    test_agent_has_asthma()
end
