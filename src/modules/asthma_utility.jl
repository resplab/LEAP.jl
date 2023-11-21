"""
    Utility

A struct for the satisfaction, or utility, that a person has given their asthma status.

For example, having asthma is worse than not having asthma, so this would have a lower utility value.
Likewise, more asthma exacerbations is worse than less asthma exacerbations, so this would
also have a lower utility value.

# Fields
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `exac`, `control`, `eq5d`.

"""
struct Utility <: Utility_Module
    parameters::Union{AbstractDict, Nothing}
end


"""
    process(agent, utility)

Compute the economic utility due to asthma.

If the person (agent) does not have asthma, the utility is the baseline utility.
If the person (agent) does have asthma, compute the reduced utility due to asthma exacerbations and
asthma control.

# Arguments
- `agent::Agent`: An agent, see [`Agent`](@ref).
- `utility::Utility`: A utility object, see [`Utility`](@ref).

"""
function process(agent::Agent, utility::Utility)
    baseline = utility.parameters[:eq5d][(agent.age, Int(agent.sex))].eq5d[1]
    if !agent.has_asthma
        return baseline
    else
        util_exacerbation = sum(agent.exac_sev_hist[1] .* utility.parameters[:exac])
        util_control = sum(agent.control .* utility.parameters[:control])
        return max(0, (baseline - util_exacerbation - util_control))
    end
end
