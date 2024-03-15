using SpecialFunctions


"""
    ExacerbationSeverity

A struct containing information about asthma exacerbation severity.
There are four levels of asthma exacerbation severity:
    1 = mild
    2 = moderate
    3 = severe
    4 = very severe

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
    function ExacerbationSeverity(config::Union{AbstractDict, Nothing})
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        hyperparameters[:α] = Array{Float64, 1}(hyperparameters[:α])
        parameters[:p] = Array{Float64, 1}(parameters[:p])
        new(hyperparameters, parameters)
    end
    function ExacerbationSeverity(
        hyperparameters::Union{AbstractDict, Nothing}, parameters::Union{AbstractDict, Nothing}
    )
        new(hyperparameters, parameters)
    end
end


"""
    ExacerbationSeverityHist

A struct containing information about the history of asthma exacerbations by severity.
There are four levels of asthma exacerbation severity:
    1 = mild
    2 = moderate
    3 = severe
    4 = very severe

# Fields
- `current_year::Array{Integer, 1}`: An array of 4 integers indicating the number of
    exacerbations for that severity level in the current year.
- `prev_year::Array{Integer, 1}`: An array of 4 integers indicating the number of
    exacerbations for that severity level in the previous year.
"""
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
- `agent::AgentModule`:: A person in the simulation, see [`Agent`](@ref).
- `exac_severity::ExacerbationSeverity`: Exacerbation severity module, see
    [`ExacerbationSeverity`](@ref).
- `control::ControlModule`: Asthma control module, see
    [`Control`](@ref).
- `exacerbation::ExacerbationModule`: Asthma exacerbation module, see
    [`Exacerbation`](@ref).

# Returns
- `Integer`: the binary probability of a hospitalization.
"""
function compute_hospitalization_prob(agent::AgentModule, exac_severity::ExacerbationSeverity,
    control::ControlModule, exacerbation::ExacerbationModule)::Integer

    max_age = agent.age - 2
    sex = agent.sex

    if max_age < 3
        return 0
    else
        # extract all the mean parameters
        # control => exac =>  sum up all the mean parameters for each tmp age
        cal_year = agent.cal_year - (agent.age - agent.asthma_age)
        total_rate = 0
        for age in agent.asthma_age:max_age
            control_levels = compute_control_levels(control, sex, age)
            # exac mean
            total_rate += compute_num_exacerbations(age, sex, cal_year, control_levels,
                exacerbation)
            cal_year += 1
        end
        # toss a coin: avg chance of having at least one hosp
        zero_prob = (
            1 / SpecialFunctions.gamma(total_rate + 1) *
            (SpecialFunctions.gamma(total_rate + 1 - exac_severity.parameters[:p][4]) /
             SpecialFunctions.gamma(1 - exac_severity.parameters[:p][4]))
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
