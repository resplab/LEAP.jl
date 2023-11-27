using Test
using LEAP
using DataFrames

@testset "test create_event_dict" begin
    simulation = LEAP.set_up(111, "CA", 2001, 40, 10, "M3");
    until_all_die = false
    max_age = simulation.max_age
    min_cal_year = simulation.starting_calendar_year
    max_cal_year = min_cal_year + simulation.time_horizon - 1
    cal_years = min_cal_year:max_cal_year
    event_dict = LEAP.create_event_dict(until_all_die, cal_years, min_cal_year,
       max_cal_year, max_age)

    @test size(event_dict["control"]) == (40, 112, 2, 3)
    @test size(event_dict["exacerbation_by_severity"]) == (40, 112, 2, 4)
    @test typeof(event_dict["asthma_prevalence_contingency_table"]) == GroupedDataFrame{DataFrame}
end

@testset "test simulation" begin
    simulation = LEAP.set_up(111, "CA", 2001, 40, 10, "M3");
    # run_test= LEAP.process(simulation,1,false,false);
end
