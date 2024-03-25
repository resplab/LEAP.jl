using Test
using LEAP
using DataFrames
using Setfield

@testset "test get_census_division_from_lat_lon" begin
    latitude = 49.262580
    longitude = -123.118720
    year = 2021
    shapefile_path = joinpath(
        LEAP.PROCESSED_DATA_PATH,
        "census_divisions/census_division_boundaries/lcd_000b21a_e.shp"
    )
    metadata_path = joinpath(
        LEAP.PROCESSED_DATA_PATH,
        "census_divisions/census_division_boundaries/lcd_000b21a_e_metadata.json"
    )
    census_boundaries = LEAP.load_census_boundaries(shapefile_path, metadata_path, year)
    census_division = LEAP.get_census_division_from_lat_lon(
        longitude=longitude,
        latitude=latitude,
        census_boundaries=census_boundaries
    )
    @test census_division.name == "Greater Vancouver"
    @test census_division.cduid == 5915
    @test census_division.year == 2021
end


@testset "test point_in_polygon" begin
    year = 2021
    shapefile_path = joinpath(
        LEAP.PROCESSED_DATA_PATH,
        "census_divisions/census_division_boundaries/lcd_000b21a_e.shp"
    )
    metadata_path = joinpath(
        LEAP.PROCESSED_DATA_PATH,
        "census_divisions/census_division_boundaries/lcd_000b21a_e_metadata.json"
    )
    census_boundaries = LEAP.load_census_boundaries(shapefile_path, metadata_path, year)
    polygon = census_boundaries.shapefile_data[260, :].geometry

    @test point_in_polygon((4017906.677490763, 2005298.410868233), polygon) == true
    @test point_in_polygon((4017906.677490763, 1005298.410868233), polygon) == false
end
