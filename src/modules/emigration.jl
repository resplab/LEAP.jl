# include("abstractModule.jl")
# include("../utils.jl")

# estimated num of newborn + proprotion of male
# birth_project = CSV.read("../processed_data/brith_projection.csv",DataFrame)

struct Emigration <: EmigrationModule
    projected_rate
    age_distribution
    table
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
