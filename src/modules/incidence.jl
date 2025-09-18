"""
    Incidence

A struct containing information about asthma incidence.

# Fields
- `hyperparameters::AbstractDict`: A dictionary containing the hyperparameters used
    to compute `β0` from a normal distribution:
    `β0_μ`: Float64, the mean of the normal distribution.
    `β0_σ`: Float64, the standard deviation of the normal distribution.
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `β0::Float64`: a constant parameter, randomly selected from a normal distribution with mean
        `β0_μ` and standard deviation `β0_σ`. See `hyperparameters`.
    `βsex::Float64`: the parameter for the sex term,
        i.e. βsex * sex
    `βage::Vector{Float64}`: an array of 5 parameters to be multiplied by functions of age,
        i.e. βage1 * f1(age) + βage2 * f2(age) + βage3 * f3(age) + βage4 * f4(age) + βage5 * f5(age)
        See `poly_age_calculator`.
    `βyear::Float64`: the parameter for the year term,
        i.e. βyear * year
    `βsexage::Vector{Float64}`: an array of 5 parameters to be multiplied by the sex and
        functions of age,
        i.e. βsexage1 * f1(age) * sex + βsexage2 * f2(age) * sex + βsexage3 * f3(age) * sex +
             βsexage4 * f4(age) * sex + βsexage5 * f5(age) * sex
        See `poly_age_calculator`.
    `βsexyear::Float64`: the parameter to be multiplied by sex and year,
        i.e. βsexyear * year * sex
    `βfam_hist::Vector{Float64}`: an array of 2 parameters to be multiplied by functions of age,
        See `log_OR_family_history`.
    `βabx_exp::Vector{Float64}`: an array of 3 parameters to be multiplied by functions of age
        and antibiotic exposure. See `log_OR_abx_exposure`.
- `min_year::Integer`: TODO.
- `max_year::Integer`: TODO.
- `max_age::Integer`: TODO.
- `correction_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
    by year, sex, and age. Each data frame contains the following columns:
        `year`: integer year.
        `sex`: sex of person, 1 = male, 0 = female.
        `age`: integer age.
        `correction`: Float64, TODO.
    See `master_occurrence_correction`.
"""
struct Incidence <: IncidenceModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    min_year::Integer
    max_year::Integer
    max_age::Integer
    correction_table::GroupedDataFrame{DataFrame}
    function Incidence(config::AbstractDict)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:βage] = Array{Float64, 1}(parameters[:βage])
        parameters[:βsexage] = Array{Float64, 1}(parameters[:βsexage])
        parameters[:βfam_hist] = Array{Float64, 1}(parameters[:βfam_hist])
        parameters[:βabx_exp] = Array{Float64, 1}(parameters[:βabx_exp])
        max_age = config["max_age"]
        correction_table = load_incidence_correction_table()
        min_year = collect(keys(correction_table)[1])[1]+1
        max_year = collect(keys(correction_table)[length(correction_table)])[1]
        new(
            hyperparameters, parameters, min_year, max_year, max_age,
            correction_table
        )
    end
    function Incidence(hyperparameters::AbstractDict, parameters::AbstractDict,
        min_year::Integer, max_year::Integer, max_age::Integer,
        correction_table::GroupedDataFrame{DataFrame}
    )
        new(
            hyperparameters, parameters, min_year, max_year, max_age,
            incidence_correction_table
        )
    end
end


function load_incidence_correction_table()
    master_occurrence_correction = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_occurrence_correction.csv"),
        DataFrame
    )
    incidence_correction_table = groupby(
        select(
            filter([:type] => (x) -> x == "inc", master_occurrence_correction),
            Not([:type])
        ),
        [:year, :sex, :age]
    )
    return incidence_correction_table
end


"""
    agent_has_asthma(agent, incidence, inc_or_prev)

Determine whether the agent obtains a new asthma diagnosis based on age and sex.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `incidence::Incidence`: Incidence module, see [`Incidence`](@ref).
- `age::Union{Integer, Nothing}`: integer age of the agent.
- `year::Union{Integer, Nothing}`: integer year.
"""
function agent_has_asthma(
    agent::Agent, incidence::Incidence, prevalence::PrevalenceModule,
    age::Union{Integer, Nothing}=nothing,
    year::Union{Integer, Nothing}=nothing
)

    if isnothing(age)
        age = min(agent.age, incidence.max_age)
    end
    if isnothing(year)
        year = agent.year
    end

    # assume no asthma if age < 3
    if age < 3
        has_asthma = false
    elseif age == 3
        has_asthma = rand(Bernoulli(prevalence_equation(
            agent.sex, age, year, agent.has_family_hist, agent.num_antibiotic_use, prevalence
        )))
    else
        has_asthma = rand(Bernoulli(incidence_equation(
            agent.sex, age, year, agent.has_family_hist, agent.num_antibiotic_use, incidence
        )))
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
function compute_asthma_age(agent::Agent, incidence::Incidence, prevalence::PrevalenceModule,
    current_age::Integer, max_asthma_age::Integer=110
)
    # obtain the previous incidence
    min_year = incidence.min_year
    max_year = incidence.max_year
    if current_age == 3
        return 3
    else
        find_asthma_age = true
        asthma_age = 3
        tmp_sex = Int(agent.sex)
        tmp_abx_num = min(agent.num_antibiotic_use, 3)
        tmp_year = min(max(agent.year - current_age + asthma_age, min_year), max_year)
        while find_asthma_age && asthma_age < max_asthma_age
            has_asthma = agent_has_asthma(agent, incidence, prevalence, asthma_age, tmp_year)
            if has_asthma
                return asthma_age
            end
            asthma_age += 1
            asthma_age = min(asthma_age, incidence.max_age)
            tmp_year += 1
            tmp_year = min(tmp_year, max_year)
        end
        return asthma_age
    end
end


function crude_incidence(sex, age, year::Integer, parameters::AbstractDict)
    poly_age = poly_age_calculator(age)
    return exp(
        parameters[:β0] +
        parameters[:βsex] * sex +
        parameters[:βyear] * year +
        parameters[:βsexyear] * sex * year +
        sum(parameters[:βage] .* poly_age) +
        sum(parameters[:βsexage] .* sex .* poly_age)
    )
end


function incidence_equation(
    sex, age::Integer, year::Integer, has_family_hist::Bool, dose::Integer, incidence::Incidence
)
    parameters = incidence.parameters
    correction_year = min(year, incidence.max_year + 1)
    year = min(year, incidence.max_year)
    p0 = crude_incidence(sex, age, year, parameters)
    p = logistic(
        logit(p0) +
        has_family_hist * log_OR_family_history(age, parameters[:βfam_hist]) +
        log_OR_abx_exposure(age, dose, parameters[:βabx_exp]) +
        incidence.correction_table[(correction_year, sex, min(age, 63))].correction[1]
    )
    return p
end


function log_OR_family_history(age::Integer, βfam_hist::Array{Float64, 1})
    βfam_hist[1] + (min(5, age) - 3) * βfam_hist[2]
end


function log_OR_abx_exposure(age::Integer, dose::Integer, βabx_exp::Array{Float64, 1})
    if (age > 7) | (dose == 0)
        return 0
    else
        return βabx_exp[1] + βabx_exp[2] * min(age, 7) + βabx_exp[3] * min(dose, 3)
    end
end


function poly_age_calculator(
    age::Integer,
    alpha::Array{Float64, 1}=[32.07692, 32.42755, 32.76123, 32.80415, 32.54075],
    nd::Array{Float64, 1}=[
        1, 520, 179636.923076923, 47536813.3328764, 11589923664.2537,
        2683688761696.54, 594554071731935
    ]
)
    fs = zeros(6)
    fs[1] =  1 / sqrt(nd[2])
    fs[2] = (age - alpha[1]) / sqrt(nd[3])
    for i in 2:5
        fs[i+1] = (
            (age-alpha[i]) * sqrt(nd[i+1]) * fs[i] -
            nd[i+1] / sqrt(nd[i]) * fs[i-1]
        ) / sqrt(nd[i+2])
    end
    popfirst!(fs)
    fs
end
