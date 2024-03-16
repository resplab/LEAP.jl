
struct Immigration <: ImmigrationModule
    sex_ratio
    estimate
    age_distribution
    overall_rate
    table::GroupedDataFrame{DataFrame}
    function Immigration(starting_year::Integer, province::String, population_growth_type::String)
        immigration_table = load_immigration_table(starting_year, province, population_growth_type)
        new(nothing, nothing, nothing, nothing, immigration_table)
    end
    function Immigration(sex_ratio, estimate, age_distribution, overall_rate, table::GroupedDataFrame)
        new(sex_ratio, estimate, age_distribution, overall_rate, table)
    end
end


function load_immigration_table(
    starting_year::Integer, province::String, population_growth_type::String
)

    master_immigration_table = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_immigration_table_modified.csv"),
        DataFrame
    )
    immigration_table = groupby(
        select(
            select(
                filter(
                    ([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year
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
    get_num_new_immigrants(immigration, num_new_born, cal_year_index)

Get the number of new immigrant to Canada in a given year.

# Arguments
- `immigration::Immigration`: an Immigration object, see  [`Immigration`](@ref).
- `num_new_born::Integer`: number of births in the given year of the simulation.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.

# Returns
- `Integer`: the number of new immigrants to Canada in a given year.
"""
function get_num_new_immigrants(immigration::Immigration, num_new_born::Integer, cal_year_index::Integer)
    num_new_immigrants = ceil(
        Int,
        num_new_born * sum(immigration.table[cal_year_index].n_prop_birth)
    )
    return num_new_immigrants
end
