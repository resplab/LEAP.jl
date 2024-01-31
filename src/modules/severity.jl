using SpecialFunctions


"""
    ExacerbationSeverity

A struct containing information about asthma exacerbation severity.
There are four levels of asthma exacerbation severity: mild, moderate, severe, and very severe.

# Fields
- `hyperparameters::Union{AbstractDict, Nothing}`: A dictionary containing the hyperparameters used
    in the Dirichlet-multinomial distribution. See
    https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.Dirichlet.
    `k`: integer, number of trials.
    `α`: parameter vector, length 4.
- `parameters::Union{AbstractDict, Nothing}`: A dictionary containing the following keys:
    `p`: a probability vector giving the probability of each exacerbation type, using the
        Dirichlet-multinomial distribution.
    `βprev_hosp_ped`: Float64, parameter for previous hospitalizations due to asthma in childhood.
    `βprev_hosp_adult`: Float64, parameter for previous hospitalizations due to asthma in adulthood.
"""
struct ExacerbationSeverity <: ExacerbationSeverityModule
    hyperparameters::Union{AbstractDict, Nothing}
    parameters::Union{AbstractDict, Nothing}
end

struct ExacerbationSeverityHist <: ExacerbationSeverityHistModule
    current_year::Array{Integer, 1}
    prev_year::Array{Integer, 1}
end


"""
    compute_distribution_exac_severity(exac_severity, num_current_year, prev_hosp, age)

Compute the exacerbation severity distribution for a patient in a given year using the
Dirichlet probability vector `p` in the Multinomial distribution. See:
https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.Multinomial.

For example, if the patient has `num_current_year` = 10 exacerbations in the current year,
then the output might be:

mild | moderate | severe | very severe
2    | 1        | 6      | 1

# Arguments
- `exac_severity::ExacerbationSeverity`: Exacerbation severity parameters, see
    [`ExacerbationSeverity`](@ref).
- `num_current_year`: the number of asthma exacerbations the patient has had this year. Will be used
    as the number of trials in the Multinomial distribution.
- `prev_hosp::Bool`: has patient been previously hospitalized for asthma?
- `age::Integer`: the age of the person in years.

# Returns
- `Array{Integer, 1}`: the distribution of asthma exacerbations by exacerbation type for the
    current year.
"""
function compute_distribution_exac_severity(exac_severity::ExacerbationSeverity,
    num_current_year::Integer, prev_hosp::Bool, age::Integer)::Array{Integer, 1}

    p = copy(exac_severity.parameters[:p])
    index_very_severe = 4
    index_max = index_very_severe

    if num_current_year == 0
        return zeros(index_max)
    else
        if prev_hosp
            weight = copy(p[1:3])
            weight = weight / sum(weight)
            p[index_very_severe] = p[index_very_severe] * (
                age < 14
                ? exac_severity.parameters[:βprev_hosp_ped]
                : exac_severity.parameters[:βprev_hosp_adult]
            )
            p[1:3] .= weight * (1 - p[index_very_severe])
        end

        return rand(Multinomial(num_current_year, p))
    end
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
    sex = sim.agent.sex

    if max_age < 3
        return 0
    else
        # extract all the mean parameters
        # control => exac =>  sum up all the mean parameters for each tmp age
        tmp_cal_year = sim.agent.cal_year - (sim.agent.age - asthma_age)
        total_rate = 0
        for tmp_age in asthma_age:max_age
            control_levels = compute_control_levels(sim.control, sex, tmp_age)
            # exac mean
            total_rate += compute_num_exacerbations(tmp_age, sex, tmp_cal_year, control_levels,
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


"""
    random_parameter_initialization!(exac_severity)

Compute the probability vector `p` from the Dirichlet distribution. See:
https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.Dirichlet.

Mutates the exac_severity object.

# Arguments
- `exac_severity::ExacerbationSeverity`: Exacerbation severity parameters, see
    [`ExacerbationSeverity`](@ref).
"""
function random_parameter_initialization!(exac_severity::ExacerbationSeverity)
    exac_severity.parameters[:p] = rand(
        Dirichlet(exac_severity.hyperparameters[:α] * exac_severity.hyperparameters[:k])
    )
end
