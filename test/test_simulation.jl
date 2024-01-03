using Test
using LEAP
using DataFrames


@testset "test simulation" begin
    simulation = LEAP.set_up(111, "CA", 2001, 40, 10, "M3");
    LEAP.process(simulation, 1, false, false);
end
