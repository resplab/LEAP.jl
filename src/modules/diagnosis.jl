"""
    Diagnosis

A struct containing information about asthma diagnosis and misdiagnosis.

# Fields
- `true_positive_rates::Union{GroupedDataFrame{DataFrame}, Nothing}`: A data frame grouped
    by year. Gives the probability of having asthma given that you have been diagnosed with
    asthma. Each data frame contains the following columns:
        `year`: integer year.
        `age`: integer age.
        `F`: Float64 - probability of having asthma given being diagnosed with asthma, female.
        `M`: Float64 - probability of having asthma given being diagnosed with asthma, male.
        `province`: A string indicating the province abbreviation, e.g. "BC".
    See `master_dx`.
- `false_negative_rates::Union{GroupedDataFrame{DataFrame}, Nothing}`: A data frame grouped
    by year. Gives the probability of having asthma given that you have not been diagnosed with
    asthma. Each data frame contains the following columns:
        `year`: integer year.
        `age`: integer age.
        `F`: Float64 - probability of having asthma given not being diagnosed with asthma, female.
        `M`: Float64 - probability of having asthma given not being diagnosed with asthma, male.
        `province`: A string indicating the province abbreviation, e.g. "BC".
    See `master_mis_dx`.
"""
struct Diagnosis <: DiagnosisModule
    true_positive_rates::Union{GroupedDataFrame{DataFrame}, Nothing}
    false_negative_rates::Union{GroupedDataFrame{DataFrame}, Nothing}
end


"""
    agent_has_asthma(agent, diagnosis)

Determine whether the agent has a correct diagnosis based on age and sex.
If the agent has asthma, adjust for false positives.
If the agent does not have asthma, adjust for false negatives.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `diagnosis::Diagnosis`: Diagnosis module, see [`Diagnosis`](@ref).
"""
function agent_has_asthma(agent::Agent, diagnosis::Diagnosis)
    max_year_index = length(diagnosis.true_positive_rates)
    age_index = agent.age - 3

    if agent.sex == 0
        column = "F"
    else
        column = "M"
    end

    if agent.age < 4
        has_asthma = agent.has_asthma
    else
        cal_year_index = min(max_year_index, agent.cal_year_index)
        if agent.has_asthma # if yes asthma, can get mis dx of not having asthma
            true_positive_prob = diagnosis.true_positive_rates[cal_year_index][age_index, column]
            has_asthma = rand(Bernoulli(true_positive_prob))
        else # if no asthma but then can get mis dx of asthma
            false_negative_prob = diagnosis.false_negative_rates[cal_year_index][age_index, column]
            has_asthma = rand(Bernoulli(false_negative_prob))
        end
    end
    return has_asthma
end
