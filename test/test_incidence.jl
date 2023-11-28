using Test
using LEAP
using DataFrames

@testset "test agent_has_asthma" begin
    simulation = LEAP.set_up(111, "CA", 2001, 40, 10, "M3");
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = process_birth(
        cal_year, tmp_cal_year_index, simulation.birth,
        simulation.antibioticExposure, simulation.familyHistory
    )
    @test agent.has_asthma == false
    agent_has_asthma(agent, simulation.incidence)

end
