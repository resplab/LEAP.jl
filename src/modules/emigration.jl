
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
                    ([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year
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
    compute_prob_emigration(agent, death)

Determine the probability of emigration of an agent (person) in a given year.

# Arguments
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.
- `age::Integer`: age of the person.
- `sex::Bool`: sex of the person, 1 = male, 0 = female.
- `emigration::Emigration`: Emigration module, see [`Emigration`](@ref).
"""
function compute_prob_emigration(cal_year_index::Integer, age::Integer, sex::Bool,
    emigration::Emigration)

    if age == 0
        return false
    else
        return rand(Bernoulli(emigration.table[cal_year_index][min(age, 100), sex+3]))
    end
end
