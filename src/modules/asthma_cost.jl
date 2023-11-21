"""
    Cost

A struct for the financial cost of asthma.

# Fields
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `exac`, `control`.

"""
struct Cost <: Cost_Module
    parameters::Union{AbstractDict, Nothing}
end


"""
    process(agent, cost)

Compute the financial cost due to asthma.

If the person (agent) does not have asthma, the cost is $0.
If the person (agent) does have asthma, compute the cost due to asthma exacerbations
and asthma control.

# Arguments
- `agent::Agent`: An agent, see [`Agent`](@ref).
- `cost::Cost`: A cost object, see [`Cost`](@ref).

"""
function process(agent::Agent, cost::Cost)
    if !agent.has_asthma
        return 0
    else
        return (
            sum(agent.exac_sev_hist[1] .* cost.parameters[:exac]) +
            sum(agent.control .* cost.parameters[:control]
        )
    end
end
