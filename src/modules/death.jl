# include("abstractmodule.jl")
# include("utils")

# prob of death by age and sex
# life_table = CSV.read("../processed_data/life_table.csv",DataFrame);

"""
    Death

A struct containing information about the probability of death for a given agent in a given year.

# Fields
- `parameters::AbstractDict`: A dictionary containing three keys,
    `:β0`, `:β1`, `:β2`.
- `life_table`: TODO.
"""
struct Death <: DeathModule
    parameters::AbstractDict
    life_table
    function Death(config::AbstractDict, province::String, starting_year::Integer)
        parameters = string_to_symbols_dict(config["parameters"])
        master_life_table = load_life_table()
        life_table = groupby(select(unstack(
            select(
                select(
                    filter((
                        [:year, :province] => (x, y) -> x >= starting_year
                        && y == province
                        ),
                        master_life_table
                    ),
                Not(:se)),
            Not(:province)),
        :sex, :prob_death), :F, :M, :year), :year)
        new(parameters, life_table)
    end
    function Death(parameters::AbstractDict, life_table)
        new(parameters, life_table)
    end
end


function load_life_table()
    master_life_table = CSV.read(joinpath(PROCESSED_DATA_PATH, "master_life_table.csv"), DataFrame)
    return master_life_table
end


"""
    compute_prob_death(agent, death)

Determine the probability of death of an agent in a given year.

# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `death::Death`: Death module, see [`Death`](@ref).
"""
function compute_prob_death(agent::Agent, death::Death)
    p = death.life_table[agent.cal_year_index][agent.age+1, agent.sex+1]
    if p == 1
        return true
    end
    # calibration
    or = p/(1-p)*exp(
        death.parameters[:β0] +
        death.parameters[:β1]*agent.cal_year_index +
        death.parameters[:β2]*agent.age
    )
    p = max(min(or/(1 + or), 1), 0)
    return rand(Bernoulli(p))
end
