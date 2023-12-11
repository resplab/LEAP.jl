using SpecialFunctions

struct ExacerbationSeverity <: Exacerbation_Severity_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
end

function process_severity(exac_severity::ExacerbationSeverity, num::Integer, prev_hosp::Bool,
    age::Integer)
    tmp_p = copy(exac_severity.parameters[:p])
    len_p = length(tmp_p)

    if isnothing(num) || num==0
        return zeros(len_p)
    else

        if prev_hosp
            tmp_weight = copy(tmp_p[1:3])
            tmp_weight = tmp_weight / sum(tmp_weight)
            tmp_p[len_p] = tmp_p[len_p]*(
                age < 14 ? exac_severity.parameters[:βprev_hosp_ped] : exac_severity.parameters[:βprev_hosp_adult])
            tmp_p[1:(len_p-1)] .= tmp_weight * (1-tmp_p[len_p])
        end

        return rand(Multinomial(num,tmp_p))
    end
    # rescale p
end

# CIHI Data on Quebec:

# 1) Impute the fraction of admissions due to asthma for Quebec (based on trend analysis) for prior to April 2001 and re-do the analyses;
# 2) Re-do the analyses excluding Quebec before 2002;
# 3) Re-do the analyses for a suitable period during which Quebec reported the primary diagnosis type (from April 2001 onwards).


function process_ctl(age,sex,ctl::Control)
    age_scaled = age / 100
    function control_prediction_here(eta::Float64,theta::Union{Float64,Vector{Float64}};inv_link::Function=StatsFuns.logistic)::Union{Float64,Vector{Float64}}
        theta = [-1e5;theta;1e5]
        [inv_link(theta[j+1] - eta) - inv_link(theta[j] - eta) for j in 1:(length(theta)-1)]
    end

    control_prediction_here(ctl.parameters[:β0]+
    age_scaled*ctl.parameters[:βage]+
    sex*ctl.parameters[:βsex]+
    age_scaled * sex * ctl.parameters[:βsexage] +
    age_scaled^2 * sex * ctl.parameters[:βsexage2] +
    age_scaled^2 * ctl.parameters[:βage2], ctl.parameters[:θ])
end


"""
    compute_hospitalization_prob(agent, incidence, current_age)

Determine whether a person (agent) has been hospitalized due to an asthma exacerbation.

https://stats.stackexchange.com/questions/174952/marginal-probability-function-of-the-dirichlet-multinomial-distribution

# Arguments
- `exac_severity::ExacerbationSeverity`: Exacerbation severity module, see  [`ExacerbationSeverity`](@ref).
- `asthma_age::Integer`: The age when the person (agent) was first diagnosed with asthma.

# Returns
- `Bool`: the binary probability of a hospitalization.
"""
function compute_hospitalization_prob(exac_severity::ExacerbationSeverity, asthma_age::Integer, sim)
    max_age = sim.agent.age - 2
    tmp_sex = sim.agent.sex

    if max_age < 3
        return 0
    else
        # extract all the mean parameters
        # control => exac =>  sum up all the mean parameters for each tmp age
        tmp_cal_year = sim.agent.cal_year - (sim.agent.age - asthma_age)
        total_rate = 0
        for tmp_age in asthma_age:max_age
            tmp_control = process_ctl(tmp_age,tmp_sex,sim.control);
            # exac mean
            total_rate += compute_num_exacerbations(tmp_age, tmp_sex, tmp_cal_year, tmp_control,
                sim.exacerbation)
            tmp_cal_year +=1
        end
        # toss a coin: avg chance of having at least one hosp
        zero_prob = (
            1 / SpecialFunctions.gamma(total_rate+1) *
            (SpecialFunctions.gamma(total_rate+1-exac_severity.parameters[:p][4]) /
             SpecialFunctions.gamma(1-exac_severity.parameters[:p][4]))
        )
        p = 1 - min(max(zero_prob, 0), 1)
        return Int(rand(Bernoulli(p)))
    end
end


function random_parameter_initialization!(exac_severity::ExacerbationSeverity)
    exac_severity.parameters[:p] = rand(Dirichlet(exac_severity.hyperparameters[:p0_μ]*exac_severity.hyperparameters[:p0_σ]))
end
