struct FamilyHistory <: FamilyHistory_Module
    hyperparameters
    parameters::AbstractDict
end

function process(ag::Agent,fam::FamilyHistory)
    rand(Bernoulli(fam.parameters[:p]))
end

function process_initial(ag::Agent,fam::FamilyHistory)
    rand(Bernoulli(fam.parameters[:p]))
end
