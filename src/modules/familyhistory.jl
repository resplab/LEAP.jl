struct FamilyHistory <: FamilyHistoryModule
    parameters::AbstractDict
    function FamilyHistory(config::AbstractDict)
        parameters = string_to_symbols_dict(config["parameters"])
        new(parameters)
    end
    function FamilyHistory(parameters::AbstractDict)
        new(parameters)
    end
end

function has_family_history_of_asthma(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end
