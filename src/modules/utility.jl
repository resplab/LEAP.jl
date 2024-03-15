"""
    Utility

A struct containing information about the disutility from having asthma.

# Fields
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `control`: A vector of numbers.
    `exac`: A vector of numbers.
    `eq5d`: TODO.
"""
struct Utility <: UtilityModule
    parameters::Union{AbstractDict, Nothing}
    function Utility(config::Union{AbstractDict, Nothing})
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:exac] = Array{Float64, 1}(parameters[:exac])
        parameters[:control] = Array{Float64, 1}(parameters[:control])
        parameters[:eq5d] = load_eq5d()
        new(parameters)
    end
    function Utility(parameters::Union{AbstractDict, Nothing})
        new(parameters)
    end
end


function load_eq5d()
    eq5d = CSV.read(joinpath(PROCESSED_DATA_PATH, "eq5d_canada.csv"), DataFrame)
    eq5d = groupby(eq5d, [:age, :sex])
    return eq5d
end


"""
    compute_cost(agent, incidence)

Compute the utility for the current year due to asthma exacerbations and control. If the agent
(person) doesn't have asthma, return the baseline utility.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `utility::Utility`: Utility module, see [`Utility`](@ref).
"""
function compute_utility(agent::Agent, utility::Utility)
    baseline = utility.parameters[:eq5d][(agent.age, Int(agent.sex))].eq5d[1]
    if !agent.has_asthma
        return baseline
    else
        disutil_exac = sum(agent.exac_sev_hist.current_year .* utility.parameters[:exac])
        disutil_control = sum(agent.control_levels[:as_array] .* utility.parameters[:control])
        return max(0, (baseline - disutil_exac - disutil_control))
    end
end
