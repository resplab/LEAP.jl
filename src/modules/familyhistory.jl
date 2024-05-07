"""
    FamilyHistory

A struct containing information about family history of asthma.

# Fields
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `p`: Float64, the probability that an agent has a family history of asthma.
"""
struct FamilyHistory <: FamilyHistoryModule
    parameters::AbstractDict
    function FamilyHistory(config::AbstractDict)
        parameters = string_to_symbols_dict(config["parameters"])
        if parameters[:p] > 1 || parameters[:p] < 0
            throw(ArgumentError(
                "p must be a probability between 0 and 1, received $(parameters[:p])."
            ))
        end
        new(parameters)
    end
    function FamilyHistory(parameters::AbstractDict)
        if parameters[:p] > 1 || parameters[:p] < 0
            throw(ArgumentError(
                "p must be a probability between 0 and 1, received $(parameters[:p])."
            ))
        end
        new(parameters)
    end
end


"""
    has_family_history_of_asthma(family_history)

Use a Bernoulli distribution to determine whether an agent has a family history of asthma.

# Arguments
- `family_history::FamilyHistory`: A family history struct, see [`FamilyHistory`](@ref).

# Returns
- `Bool`: Whether or not an agent has a family history of asthma.
"""
function has_family_history_of_asthma(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end
