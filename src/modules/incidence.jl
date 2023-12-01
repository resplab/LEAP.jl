struct Incidence <: Incidence_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    incidence_table
    prevalence_table
    calibration_table
    min_year
    max_year
    initial_distribution
end


function agent_has_asthma(agent::Agent, incidence::Incidence, abx)
    tmp_age = min(agent.age, 95)
    # max_year = length(inc.incidence_table)
    tmp_year = min(agent.cal_year, incidence.max_year)
    if tmp_age < 3
        has_asthma = false
    elseif tmp_age < 7
        try
            has_asthma = rand(Bernoulli(
                incidence.calibration_table[(
                    tmp_year, Int(agent.sex), Int(agent.family_hist), min(agent.num_antibiotic_use,3)
                )][tmp_age - 2, "calibrated_inc"]))
        catch
            println(tmp_year, " ", agent.sex, " ", agent.family_hist, " ",
                agent.num_antibiotic_use, " ", tmp_age)
        end
    else # no effect of Abx beyond 7 years of age
        rand(Bernoulli(incidence.calibration_table[(
            tmp_year, Int(agent.sex), Int(agent.family_hist), 0
        )][tmp_age - 2, "calibrated_inc"]))
    end
    return has_asthma
end


function agent_has_asthma(agent::Agent, incidence::Incidence)
    tmp_age = min(agent.age, 95) - 1
    max_year = incidence.max_year
    tmp_year = min(agent.cal_year, max_year) - 1
    # assume no asthma if age < 3
    if tmp_age < 3
        has_asthma = false
    elseif tmp_age < 7
        has_asthma = rand(Bernoulli(incidence.calibration_table[(
            tmp_year, Int(agent.sex), Int(agent.family_hist), min(agent.num_antibiotic_use, 3)
            )][tmp_age - 2, "calibrated_prev"]))
    else # no effect of Abx beyond 7 years of age
        has_asthma = rand(Bernoulli(incidence.calibration_table[(
            tmp_year, Int(agent.sex), Int(agent.family_hist), 0
            )][tmp_age - 2, "calibrated_prev"]))
    end
    return has_asthma
end

# initialization means prevalence ! ! !


"""
    compute_asthma_age(agent, incidence, current_age)

Compute the age at which the person (agent) is first diagnosed with asthma.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `incidence::Incidence`: Asthma incidence, see [`Incidence`](@ref).
"""
function compute_asthma_age(agent::Agent, incidence::Incidence, current_age::Integer)
    # obtain the previous incidence
    min_year = incidence.min_year
    max_year = incidence.max_year
    if current_age == 3
        return 3
    else
        find_asthma_age = true
        asthma_age = 3
        tmp_family_hist = Int(agent.family_hist)
        tmp_sex = Int(agent.sex)
        tmp_abx_num = min(agent.num_antibiotic_use, 3)
        tmp_year = min(max(agent.cal_year-current_age+asthma_age, min_year), max_year)
        while find_asthma_age && asthma_age < 110
            if rand(Bernoulli(
                incidence.calibration_table[(tmp_year,tmp_sex,tmp_family_hist,(asthma_age < 7 ? min(tmp_abx_num,3) : 0))][asthma_age-2,"calibrated_inc"]))
                return asthma_age
            end
            asthma_age += 1
            asthma_age = min(asthma_age,95)
            tmp_year += 1
            tmp_year = min(tmp_year,max_year)
        end
        return asthma_age
    end
end

function prevalence_logit_prob(p::Real,age::Int,param,family_hist,num_abx)
    return p
end

function incidence_logit_prob(p::Real,sex::Bool,age::Int64,cal_year::Int64,ABX::Int64,parameters::AbstractDict,abxlogOR,calibration)
    return p
end


function random_parameter_initialization!(inc::Incidence)
    inc.parameters[:β0] = rand(Normal(inc.hyperparameters[:β0_μ], inc.hyperparameters[:β0_σ]))
end
