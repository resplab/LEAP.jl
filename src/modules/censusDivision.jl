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
struct CensusDivision <: CensusDivisionModule
    federal_census_division::Union{String, Nothing}
    year::Union{Integer, Nothing}
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
