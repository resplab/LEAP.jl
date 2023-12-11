struct Reassessment <: Reassessment_Module
    table
end

function agent_has_asthma(agent::Agent, reassessment::Reassessment)
    max_year = length(reassessment.table)
    if agent.age < 4
        agent.has_asthma
    else
        rand(Bernoulli(reassessment.table[min(agent.cal_year_index,max_year)][agent.age-3, agent.sex+3]))
    end
end
