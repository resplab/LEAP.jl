using Test
using LEAP
using DataFrames
using Setfield

@testset "test load_pollution_data" begin
    pollution_data = LEAP.load_pollution_data()
    grouped_keys = keys(pollution_data)
    @test length(grouped_keys) == 4
    @test grouped_keys[1][:SSP] == "SSP1_2.6"
    @test grouped_keys[2][:SSP] == "SSP2_4.5"
    @test grouped_keys[3][:SSP] == "SSP3_7.0"
    @test grouped_keys[4][:SSP] == "SSP5_8.5"
    @test pollution_data[("SSP1_2.6",)][1, :month] == 1
    @test pollution_data[("SSP1_2.6",)][10, :CDUID] == 1001
end

@testset "test PollutionTable constructor no data" begin
    pollution_table = LEAP.PollutionTable()
    pollution_data = pollution_table.data
    @test typeof(pollution_data) == GroupedDataFrame{DataFrame}
    grouped_keys = keys(pollution_data)
    @test length(grouped_keys) == 4
    @test grouped_keys[1][:SSP] == "SSP1_2.6"
    @test grouped_keys[2][:SSP] == "SSP2_4.5"
    @test grouped_keys[3][:SSP] == "SSP3_7.0"
    @test grouped_keys[4][:SSP] == "SSP5_8.5"
    @test pollution_data[("SSP1_2.6",)][1, :month] == 1
    @test pollution_data[("SSP1_2.6",)][10, :CDUID] == 1001
end

@testset "test PollutionTable constructor data" begin
    pollution_data = LEAP.load_pollution_data()
    pollution_table = LEAP.PollutionTable(pollution_data)
    pollution_data = pollution_table.data
    @test typeof(pollution_data) == GroupedDataFrame{DataFrame}
    grouped_keys = keys(pollution_data)
    @test length(grouped_keys) == 4
    @test grouped_keys[1][:SSP] == "SSP1_2.6"
    @test grouped_keys[2][:SSP] == "SSP2_4.5"
    @test grouped_keys[3][:SSP] == "SSP3_7.0"
    @test grouped_keys[4][:SSP] == "SSP5_8.5"
    @test pollution_data[("SSP1_2.6",)][1, :month] == 1
    @test pollution_data[("SSP1_2.6",)][10, :CDUID] == 1001
end


@testset "test Pollution constructor" begin
    cduid = 5915
    year = 2028
    month = 2
    SSP = "SSP2_4.5"
    pollution = LEAP.Pollution(cduid, year, month, SSP)
    @test pollution.cduid == cduid
    @test pollution.year == year
    @test pollution.month == month
    @test pollution.SSP == SSP
    @test pollution.wildfire_pm25_scaled == 0.09851602637
    @test pollution.total_pm25 == 3.168516026
end
