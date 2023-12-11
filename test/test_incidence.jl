using Test
using LEAP
using DataFrames

@testset "test agent_has_asthma" begin
    max_age = 111
    province = "CA"
    starting_year = 2001
    time_horizon = 40
    n = 10
    population_growth_type = "M3"
    incidence = LEAP.set_up_incidence(starting_year, province)
    birth = LEAP.set_up_birth(starting_year, population_growth_type, province)
    antibiotic_exposure = LEAP.set_up_antibiotic_exposure()
    family_history = LEAP.set_up_family_history()
    cal_year = 2002
    tmp_cal_year_index = cal_year - 2001 + 1
    agent = process_birth(
        cal_year, tmp_cal_year_index, birth,
        antibiotic_exposure, family_history
    )
    @test agent.has_asthma == false
    agent_has_asthma(agent, incidence)

end
