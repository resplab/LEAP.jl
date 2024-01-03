import Base.@kwdef

"""
    CensusTable

A struct containing information about Canadian census divisions.

Please see: Statistics Canada. Table 98-10-0010-02 
Population and dwelling counts: Canada, provinces and territories, and federal electoral districts
(2013 Representation Order)"
https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=9810001002

# Fields
- `data::Union{DataFrame, Nothing}`: A data frame with the following columns:
    `federal_census_division`: the federal census division name.
    `province`: the two-letter province identifier.
    `population`: the number of residents living in the census division.
    `area_km2`: the area of the census division in kilometres squared.
    `population_density_per_square_km`: the population density per square kilometre.
- `year::Integer`: the year the census population data was collected.
"""
struct CensusTable <: CensusTableModule
    data::Union{GroupedDataFrame{DataFrame}, Nothing}
    year::Union{Integer, Nothing}
end


"""
    CensusBoundaries

A struct containing information about Canadian federal electoral district boundaries.

Please see: Statistics Canada. 2021 Census – Boundary file.
Type: Cartographic Boundary File (CBF)
Administrative Boundaries: Federal Electoral Districts, 2013 representation order
Downloadable: Shapefile (.shp)

https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21

The shapefile given is in the Lambert Conformal Conic projection. More details here:
https://www150.statcan.gc.ca/n1/en/pub/92-160-g/92-160-g2021002-eng.pdf?st=0cwqeEP1

# Fields
- `shapefile_data::Union{DataFrame, Nothing}`: A data frame with the following columns:
    `geometry`: the Shapefile.Polygon of the census boundary.
    `FEDUID`: TODO.
    `DGUID`: TODO.
    `FEDNAME`: the federal electoral district name, in French and English.
    `FEDENAME`: the federal electoral district name, in English.
    `FEDFNAME`: the federal electoral district name, in French.
    `LANDAREA`: the area of the electoral district in square kilometres.
    `PRUID`: integer, province id.
- `year::Integer`: the year the census population data was collected.
- `reference_longitude::Number`: the reference longitude.
- `reference_latitude::Number`: the reference latitude.
- `first_standard_parallel::Number`: the first standard parallel in degrees.
- `second_standard_parallel::Number`: the second standard parallel in degrees.
- `false_easting::Number`: false easting.
- `false_northing::Number`: false northing.

"""
@kwdef struct CensusBoundaries <: CensusBoundariesModule
    shapefile_data::Union{DataFrame, Nothing}
    year::Union{Integer, Nothing}
    reference_longitude::Number
    reference_latitude::Number
    first_standard_parallel::Number
    second_standard_parallel::Number
    false_easting::Number
    false_northing::Number
end


"""
    CensusDivison

A struct containing information about a Canadian census division.

Please see: Statistics Canada. Table 98-10-0010-02 
Population and dwelling counts: Canada, provinces and territories, and federal electoral districts
(2013 Representation Order)"
https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=9810001002

# Fields
- `federal_census_division::Union{String, Nothing}`: the federal census division name.
- `year::Integer`: the year the census population data was collected.
"""
@kwdef struct CensusDivision <: CensusDivisionModule
    federal_census_division::Union{String, Nothing}
    year::Union{Integer, Nothing}
end


function load_census_boundaries(shapefile_path::String, metadata_path::String, year::Integer=2021)
    table = Shapefile.Table(shapefile_path)
    metadata = JSON.parsefile(metadata_path)
    census_boundaries = CensusBoundaries(
        shapefile_data=DataFrame(table),
        year=metadata["year"],
        reference_longitude=metadata["reference_longitude"],
        reference_latitude=metadata["reference_latitude"],
        first_standard_parallel=metadata["first_standard_parallel"],
        second_standard_parallel=metadata["second_standard_parallel"],
        false_easting=metadata["false_easting"],
        false_northing=metadata["false_northing"]
    )
    return census_boundaries
end


"""
    assign_census_division(census_table, province, year)

Randomly assign a census division based on population.

# Arguments
- `census_table::CensusTable`: Struct containing Canadian census data, see  [`CensusTable`](@ref).
- `province::String`: the two-letter province identifier.
- `year::Integer`: the year the census was taken (default = 2021).

# Returns
- `CensusDivision`: the census division assigned.
"""
function assign_census_division(census_table::CensusTable, province::String, year::Integer=2021)
    if province == "CA"
        census_table_all = DataFrames.DataFrame(census_table.data)
        census_division_name = sample(
            census_table_all[!, "federal_census_division"],
            Weights(census_table_all[!, "population"])
        )
    else
        census_table_province = census_table.data[(province,)]
        census_division_name = sample(
            census_table_province[!, "federal_census_division"],
            Weights(census_table_province[!, "population"])
        )
    end
    census_division = CensusDivision(
        census_division_name,
        year
    )
    return census_division
end


"""
    get_census_division_from_lat_lon(longitude, latitude, census_boundaries)

Given a latitude and longitude, find the corresponding census division.

# Arguments

- `longitude::Float64`: the longitude.
- `latitude::Float64`: the latitude.
- `census_boundaries::CensusBoundaries`: a struct containing information about Canadian census divisions.

"""
function get_census_division_from_lat_lon(; longitude::Float64, latitude::Float64,
    census_boundaries::CensusBoundaries)::CensusDivision

    point = get_lambert_conformal_from_lat_lon(
        λ=longitude,
        ϕ=latitude,
        λ_0=census_boundaries.reference_longitude,
        ϕ_0=census_boundaries.reference_latitude,
        ϕ_1=census_boundaries.first_standard_parallel,
        ϕ_2=census_boundaries.second_standard_parallel,
        x_0=census_boundaries.false_easting,
        y_0=census_boundaries.false_northing
    )

    census_division = CensusDivision(
        federal_census_division=nothing,
        year=census_boundaries.year
    )
    is_point_in_polygon = false

    for row_index in 1:size(census_boundaries.shapefile_data)[1]
        row = census_boundaries.shapefile_data[row_index, :]
        polygon = row.geometry
        println(row_index)
        if point_in_polygon(point, polygon)
            @set! census_division.federal_census_division = row.FEDENAME
            is_point_in_polygon = true
            break
        end
    end

    if !is_point_in_polygon
        throw(error("Could not find point in any of the federal electoral districts."))
    else
        return census_division
    end
end


"""
    point_in_polygon(point, polygon)

Determine whether or not a given point is within a polygon.

# Arguments

- `point::Tuple{Number, Number}`: the point we want to check.
- `polygon::Shapefile.Polygon`: a polygon object containing the points of its boundary.

"""
function point_in_polygon(point::Tuple{Number, Number}, polygon::Shapefile.Polygon)::Bool
    polygon_array = [];
    polygon_array = GeoInterface.coordinates(polygon)
    return point_in_polygon_array(point, polygon_array)
end


function point_in_polygon_array(point::Tuple{Number, Number}, polygon_array::Array)::Bool
    first_entry = polygon_array[1]
    if typeof(first_entry) == Vector{Float64} && typeof(first_entry[1]) == Float64 && typeof(first_entry[2]) == Float64
        return PolygonOps.inpolygon(point, polygon_array)
    else
        is_in_polygon_array = false
        for polygon_sub_array in polygon_array
            if point_in_polygon_array(point, polygon_sub_array)
                is_in_polygon_array = true
                break
            end
        end
        return is_in_polygon_array
    end
end


"""
    get_lambert_conformal_from_lat_lon(λ, ϕ, λ_0, ϕ_0, ϕ_1, ϕ_2, x_0, y_0)

Given a latitude and longitude, find the Lamber Conformal Conic projection coordinates.

See: https://www.linz.govt.nz/guidance/geodetic-system/understanding-coordinate-conversions/projection-conversions/lambert-conformal-conic-geographic-transformation-formulae
also: https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980

# Arguments

- `λ::Float64`: the longitude.
- `ϕ::Float64`: the latitude.
- `λ_0::Float64`: the reference longitude.
- `ϕ_0::Float64`: the reference latitude.
- `ϕ_1::Float64`: the first standard parallel in degrees.
- `ϕ_2::Float64`: the second standard parallel in degrees.
- `x_0::Number`: false easting.
- `y_0::Number`: false northing.
"""
function get_lambert_conformal_from_lat_lon(; λ::Number, ϕ::Number, λ_0::Number, ϕ_0::Number,
    ϕ_1::Number, ϕ_2::Number, x_0::Number=0, y_0::Number=0)::Tuple{Number, Number}

    λ = deg2rad(λ)
    ϕ = deg2rad(ϕ)
    λ_0 = deg2rad(λ_0)
    ϕ_0 = deg2rad(ϕ_0)
    ϕ_1 = deg2rad(ϕ_1)
    ϕ_2 = deg2rad(ϕ_2)

    R = 6378137 # Radius of Earth

    f = 0.003352810681183637418 # flattening
    e = sqrt(2*f - f^2) # eccentricity

    m_1 = cos(ϕ_1)/(sqrt(1 - e^2*sin(ϕ_1)^2))
    m_2 = cos(ϕ_2)/(sqrt(1 - e^2*sin(ϕ_2)^2))
    t = tan(π/4 - ϕ/2) / ((1 - e*sin(ϕ))/(1 + e*sin(ϕ)))^(e/2)
    t_0 = tan(π/4 - ϕ_0/2) / ((1 - e*sin(ϕ_0))/(1 + e*sin(ϕ_0)))^(e/2)
    t_1 = tan(π/4 - ϕ_1/2) / ((1 - e*sin(ϕ_1))/(1 + e*sin(ϕ_1)))^(e/2)
    t_2 = tan(π/4 - ϕ_2/2) / ((1 - e*sin(ϕ_2))/(1 + e*sin(ϕ_2)))^(e/2)

    if ϕ_1 == ϕ_2
        n = sin(ϕ_1)
    else
        n = (log(m_1) - log(m_2))/(log(t_1) - log(t_2))
    end

    F = m_1 / (n*t_1^n)
    ρ_0 = R*F*t_0^n
    ρ = R*F*t^n

    x = x_0 + ρ*sin(n*(λ - λ_0))
    y = y_0 + ρ_0 - ρ*cos(n*(λ - λ_0))

    return x, y

end
