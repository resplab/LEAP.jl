struct FamilyHistory <: FamilyHistoryModule
    hyperparameters
    parameters::AbstractDict
end

function process_family_history(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end

function process_family_history_initial(family_history::FamilyHistory)
    return rand(Bernoulli(family_history.parameters[:p]))
end
