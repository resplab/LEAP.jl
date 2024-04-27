"""
    Emigration

A struct containing information about emigration from Canada.

# Fields
- `table::GroupedDataFrame{DataFrame}`: A dataframe grouped by year, giving the probability of
    emigration for a given age, province, sex, and growth scenario:
        `year`: integer year the range 2001 - 2065.
        `age`: integer age.
        `M`: the probability of a male emigrating.
        `F`: the probability of a female emigrating.
    See `master_emigration_table.csv`.
"""
struct Emigration <: EmigrationModule
    table::GroupedDataFrame{DataFrame}
    function Emigration(starting_year::Integer, province::String, population_growth_type::String)
        emigration_table = load_emigration_table(starting_year, province, population_growth_type)
        new(emigration_table)
    end
    function Emigration(table)
        new(table)
    end
end


"""
    load_emigration_table(starting_year, province, population_growth_type)

Load the data from the `master_emigration_table.csv`.

# Arguments
- `starting_year::Integer`: the year for the data to start at. Must be between 2001-2065.
- `province::String`: a string indicating the province abbreviation, e.g. "BC". For all of Canada,
    set province to "CA".
- `population_growth_type`: Population growth type, one of:
    ["past", "LG", "HG", "M1", "M2", "M3", "M4", "M5", "M6", FA", "SA"].
    See [Stats Canada](https://www150.statcan.gc.ca/n1/pub/91-520-x/91-520-x2022001-eng.htm).


# Returns
- `GroupedDataFrame{DataFrame}`: A dataframe grouped by year, giving the probability of
    emigration for a given age, province, sex, and growth scenario.
"""
function load_emigration_table(
    starting_year::Integer, province::String, population_growth_type::String
)::GroupedDataFrame{DataFrame}

    master_emigration_table = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "migration/master_emigration_table.csv"),
        DataFrame
    )
    emigration_table = groupby(
        select(
            select(
                filter(
                    ([:year, :province, :proj_scenario] => (x, y, z) ->
                    x > starting_year
                    && y == province
                    && z == population_growth_type),
                    master_emigration_table
                ),
                Not(:province)
            ),
            Not(:proj_scenario)
        ),
        :year
    )
    return emigration_table
end


"""
    compute_prob_emigration(cal_year_index, age, sex, emigration)

Determine the probability of emigration of an agent (person) in a given year.

# Arguments
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.
- `age::Integer`: age of the person.
- `sex::Bool`: sex of the person, 1 = male, 0 = female.
- `emigration::Emigration`: Emigration module, see [`Emigration`](@ref).
"""
function compute_prob_emigration(cal_year_index::Integer, age::Integer, sex::Bool,
    emigration::Emigration
)

    if age == 0
        return false
    else
        return rand(Bernoulli(emigration.table[cal_year_index][min(age, 100), sex + 3]))
    end
end
