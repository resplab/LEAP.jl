# include("abstractModule.jl")
# include("../utils.jl")

# estimated num of newborn + proprotion of male
# birth_project = CSV.read("../processed_data/brith_projection.csv",DataFrame)

struct Immigration <: ImmigrationModule
    sex_ratio
    estimate
    age_distribution
    overall_rate
    table
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
