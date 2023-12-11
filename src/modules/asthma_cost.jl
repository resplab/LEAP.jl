"""
    AsthmaCost

A struct containing information about the cost of asthma.

# Fields
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `control`: A vector of numbers.
    `exac`: A vector of numbers.

"""
struct AsthmaCost <: Cost_Module
    parameters::Union{AbstractDict, Nothing}
end


"""
    compute_cost(agent, incidence)

Compute the cost in dollars (CAD) for the current year due to asthma exacerbations and control.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `asthma_cost::AsthmaCost`: AsthmaCost module, see [`AsthmaCost`](@ref).
"""
function compute_cost(agent::Agent, asthma_cost::AsthmaCost)
    if !agent.has_asthma
        return 0
    else
        return (
            sum(agent.exac_sev_hist.current_year .* asthma_cost.parameters[:exac]) +
            sum(agent.control .* asthma_cost.parameters[:control])
        )
    end
end
