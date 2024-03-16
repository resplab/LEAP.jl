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

function process_family_history(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end

function process_family_history_initial(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end
