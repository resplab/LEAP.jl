struct Diagnosis <: Diagnosis_Module
    table
    table_mis
end

function agent_has_asthma(agent::Agent, diagnosis::Diagnosis)
    max_year = length(diagnosis.table)
    if agent.age < 4
        agent.has_asthma
    else
        if agent.has_asthma # if yes asthma, can get mis dx of not having asthma
            rand(Bernoulli(diagnosis.table[min(max_year, agent.cal_year_index)][agent.age - 3, agent.sex + 3]))
        else # if no asthma but then can get mis dx of asthma
            rand(Bernoulli(diagnosis.table_mis[min(agent.cal_year_index, max_year)][agent.age-3, agent.sex+3]))
        end
    end
end
