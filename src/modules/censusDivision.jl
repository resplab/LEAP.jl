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
