"""
    Immigration

A struct containing information about immigration to Canada.

# Fields
- `table::GroupedDataFrame{DataFrame}`: A dataframe grouped by year, giving the probability of
    immigration for a given age, province, sex, and growth scenario:
        `year`: integer year the range 2001 - 2065.
        `age`: integer age.
        `sex`: integer, 0 = female, 1 = male.
        `prop_immigrants_birth`: the number of immigrants relative to the number of births
            in that year. To compute the number of immigrants in a given year, multiply the number
            of births by `prop_immigrants_birth`.
        `prop_immigrants_year`: the proportion of immigrants for a given age and sex relative to
            the total number of immigrants for a given year and projection scenario.
    See `master_immigration_table.csv`.
"""
struct Immigration <: ImmigrationModule
    table::GroupedDataFrame{DataFrame}
    function Immigration(
        starting_year::Integer, province::String, population_growth_type::String, max_age::Integer
    )
        immigration_table = load_immigration_table(
            starting_year, province, population_growth_type, max_age
        )
        new(immigration_table)
    end
    function Immigration(table::GroupedDataFrame)
        new(table)
    end
end


"""
    load_immigration_table(starting_year, province, population_growth_type)

Load the data from the `master_immigration_table.csv`.

# Arguments
- `starting_year::Integer`: the year for the data to start at. Must be between 2001-2065.
- `province::String`: a string indicating the province abbreviation, e.g. "BC". For all of Canada,
    set province to "CA".
- `population_growth_type`: Population growth type, one of:
    ["past", "LG", "HG", "M1", "M2", "M3", "M4", "M5", "M6", FA", "SA"].
    See [Stats Canada](https://www150.statcan.gc.ca/n1/pub/91-520-x/91-520-x2022001-eng.htm).


# Returns
- `GroupedDataFrame{DataFrame}`: A dataframe grouped by year, giving the probability of
    immigration for a given age, province, sex, and growth scenario.
"""
function load_immigration_table(
    starting_year::Integer, province::String, population_growth_type::String, max_age::Integer=111
)

    master_immigration_table = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "migration/master_immigration_table.csv"),
        DataFrame
    )
    immigration_table = groupby(
        select(
            select(
                filter(
                    ([:age, :year, :province, :proj_scenario] => (w, x, y, z) ->
                    w <= max_age
                    && x >= starting_year
                    && y == province
                    && z == population_growth_type),
                    master_immigration_table
                ),
                Not(:province)
            ),
            Not(:proj_scenario)
        ),
        :year
    )
    return immigration_table
end


"""
    get_num_new_immigrants(immigration, num_new_born, year_index)

Get the number of new immigrant to Canada in a given year.

# Arguments
- `immigration::Immigration`: an Immigration object, see  [`Immigration`](@ref).
- `num_new_born::Integer`: number of births in the given year of the simulation.
- `year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `year_index` for 2023 is 1, for 2024 is 2, etc.

# Returns
- `Integer`: the number of new immigrants to Canada in a given year.
"""
function get_num_new_immigrants(immigration::Immigration, num_new_born::Integer, year_index::Integer)
    num_new_immigrants = ceil(
        Int,
        num_new_born * sum(immigration.table[year_index].prop_immigrants_birth)
    )
    return num_new_immigrants
end
