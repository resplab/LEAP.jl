"""
    AsthmaCost

A struct containing information about the cost of asthma.

# Fields
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `control`: A vector of numbers.
    `exac`: A vector of numbers.

"""
struct AsthmaCost <: CostModule
    parameters::AbstractDict
    function AsthmaCost(config::AbstractDict)
        parameters = string_to_symbols_dict(config["parameters"])
        exchange_rate_usd_cad = config["exchange_rate_usd_cad"]
        parameters[:exac] = Array{Float64, 1}(parameters[:exac]) * exchange_rate_usd_cad
        parameters[:control] = Array{Float64, 1}(parameters[:control]) * exchange_rate_usd_cad
        new(parameters)
    end
    function AsthmaCost(parameters::AbstractDict)
        new(parameters)
    end
end


"""
    compute_cost(agent, asthma_cost)

Compute the cost in dollars (CAD) for the current year due to asthma exacerbations and control.

If `control_levels` are not present, will default to equal probability of fully controlled,
partially controlled, and uncontrolled asthma.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `asthma_cost::AsthmaCost`: AsthmaCost module, see [`AsthmaCost`](@ref).
"""
function compute_cost(agent::Agent, asthma_cost::AsthmaCost)::Float64
    if !agent.has_asthma
        return 0.0
    else
        control_levels = agent.control_levels
        if isnothing(control_levels)
            control_levels = Dict(
                :fully_controlled => 0.3333,
                :partially_controlled => 0.3333,
                :uncontrolled => 0.3333,
                :as_array => [0.3333, 0.3333, 0.3333]
            )
            @warn "Control levels have not been set, default to $control_levels"
        end
        return (
            sum(agent.exac_sev_hist.current_year .* asthma_cost.parameters[:exac]) +
            sum(control_levels[:as_array] .* asthma_cost.parameters[:control])
        )
    end
end
