"""
    Death

A struct containing information about the probability of death for a given agent in a given year.

# Fields
- `parameters::AbstractDict`: A dictionary containing three keys,
    `:β0`, `:β1`, `:β2`.
- `life_table::GroupedDataFrame{DataFrame}`: TODO.
"""
struct Death <: DeathModule
    parameters::AbstractDict
    life_table::GroupedDataFrame{DataFrame}
    function Death(config::AbstractDict, province::String, starting_year::Integer)
        parameters = string_to_symbols_dict(config["parameters"])
        master_life_table = load_life_table()
        life_table = groupby(select(unstack(
            select(
                select(
                    filter((
                        [:year, :province] => (x, y) ->
                        x >= starting_year
                        && y == province
                        ),
                        master_life_table
                    ),
                Not(:se)),
            Not(:province)),
        :sex, :prob_death), :F, :M, :year), :year)
        new(parameters, life_table)
    end
    function Death(parameters::AbstractDict, life_table::GroupedDataFrame{DataFrame})
        new(parameters, life_table)
    end
end


function load_life_table()::DataFrame
    master_life_table = CSV.read(joinpath(PROCESSED_DATA_PATH, "master_life_table.csv"), DataFrame)
    return master_life_table
end


"""
    agent_dies(agent, death)

Determine whether or not the agent dies in a given year, based on age and sex.

# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `death::Death`: Death module, see [`Death`](@ref).
"""
function agent_dies(agent::Agent, death::Death)::Bool
    is_dead = false
    p = death.life_table[agent.year_index][agent.age + 1, agent.sex + 1]
    if p == 1
        is_dead = true
    else
        # calibration
        or = p/(1 - p) * exp(
            death.parameters[:β0] +
            death.parameters[:β1] * agent.year_index +
            death.parameters[:β2] * agent.age
        )
        p = max(min(or/(1 + or), 1), 0)
        is_dead = rand(Bernoulli(p))
    end
    return is_dead
end
