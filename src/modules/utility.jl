"""
    Utility

A struct containing information about the disutility from having asthma.

# Fields
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `eq5d`: TODO.
    `βcontrol`: A vector of 3 parameters to be multiplied by the control levels, i.e.
        βcontrol1 * fully_controlled + βcontrol2 * partially_controlled + βcontrol3 * uncontrolled
    `βexac_sev_hist`: A vector of 4 parameters to be multiplied by the exacerbation
        severity history, i.e.
        βexac_sev_hist1 * mild + βexac_sev_hist2 * moderate +
        βexac_sev_hist3 * severe + βexac_sev_hist4 * very_severe
"""
struct Utility <: UtilityModule
    parameters::AbstractDict
    function Utility(config::AbstractDict)
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:eq5d] = load_eq5d()
        new(parameters)
        parameters[:βexac_sev_hist] = Array{Float64, 1}(parameters[:βexac_sev_hist])
        parameters[:βcontrol] = Array{Float64, 1}(parameters[:βcontrol])
    end
    function Utility(parameters::AbstractDict)
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
        disutility_exac = sum(
            agent.exac_sev_hist.current_year .* utility.parameters[:βexac_sev_hist]
        )
        disutility_control = sum(agent.control_levels[:as_array] .* utility.parameters[:βcontrol])
        return max(0, (baseline - disutility_exac - disutility_control))
    end
end
