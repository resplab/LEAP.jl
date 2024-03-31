"""
    Incidence

A struct containing information about asthma incidence.

# Fields
- `hyperparameters::AbstractDict`: A dictionary containing two keys,
    `β0_μ` and `β0_σ`.
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `βage`, `βage2`,`βage3`, `βage4`, `βage5`, `βageM`, `βage2M`, `βage3M`, `βage4M`, `βage5M`,
    `βsex`, `matrix_Abx_risk`, `βcal_year`, `βcal_yearM`, `βCABE`, `β0`, `β0_correction`,
    and `β0_overall_correction`.
- `parameters_prev::AbstractDict`: TODO.
- `min_year::Integer`: TODO.
- `max_year::Integer`: TODO.
- `max_age::Integer`: TODO.
- `incidence_correction_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
    by year, sex, and age. Each data frame contains the following columns:
        `year`: integer year.
        `sex`: sex of person, 1 = male, 0 = female.
        `age`: integer age.
        `correction`: Float64, TODO.
    See `master_occurrence_correction`.
- `prevalence_correction_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
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
    parameters_prev::AbstractDict
    min_year::Integer
    max_year::Integer
    max_age::Integer
    incidence_correction_table::GroupedDataFrame{DataFrame}
    prevalence_correction_table::GroupedDataFrame{DataFrame}
    function Incidence(config::AbstractDict)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        parameters_prev = string_to_symbols_dict(config["parameters_prev"])
        parameters[:βfam_hist] = Array{Float64, 1}(parameters[:βfam_hist])
        parameters[:βabx_exp] = Array{Float64, 1}(parameters[:βabx_exp])
        parameters_prev[:βfam_hist] = Array{Float64, 1}(parameters_prev[:βfam_hist])
        parameters_prev[:βabx_exp] = Array{Float64, 1}(parameters_prev[:βabx_exp])
        max_age = config["max_age"]
        incidence_correction_table = load_incidence_correction_table()
        prevalence_correction_table = load_prevalence_correction_table()
        min_year = collect(keys(incidence_correction_table)[1])[1]+1
        max_year = collect(keys(incidence_correction_table)[length(incidence_correction_table)])[1]
        new(
            hyperparameters, parameters, parameters_prev, min_year, max_year, max_age,
            incidence_correction_table, prevalence_correction_table
        )
    end
    function Incidence(hyperparameters::AbstractDict, parameters::AbstractDict,
        parameters_prev::AbstractDict, min_year::Integer, max_year::Integer, max_age::Integer,
        incidence_correction_table::GroupedDataFrame{DataFrame},
        prevalence_correction_table::GroupedDataFrame{DataFrame}
    )
        new(
            hyperparameters, parameters, parameters_prev, min_year, max_year, max_age,
            incidence_correction_table, prevalence_correction_table
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


function load_prevalence_correction_table()
    master_occurrence_correction = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_occurrence_correction.csv"),
        DataFrame
    )
    prevalence_correction_table = groupby(
        select(
            filter([:type] => (x) -> x == "prev", master_occurrence_correction),
            Not([:type])
        ),
        [:year, :sex, :age]
    )
    return prevalence_correction_table
end


"""
    agent_has_asthma(agent, incidence, inc_or_prev)

Determine whether the agent obtains a new asthma diagnosis based on age and sex.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `incidence::Incidence`: Incidence module, see [`Incidence`](@ref).
- `inc_or_prev::String`: One of "prevalence" or "incidence"
"""
function agent_has_asthma(agent::Agent, incidence::Incidence, inc_or_prev::String)
    tmp_age = min(agent.age, 95) - 1
    max_year = incidence.max_year
    tmp_year = min(agent.cal_year, max_year) - 1

    if inc_or_prev == "prevalence"
        has_asthma = agent_has_asthma_prevalence(agent, incidence)
    elseif inc_or_prev == "incidence"
        has_asthma = agent_has_asthma_incidence(agent, incidence)
    else
        throw(ArgumentError(
            "inc_or_prev must be either 'prevalence' or 'incidence', received $inc_or_prev."
        ))
    end
    return has_asthma
end


function agent_has_asthma_incidence(
    agent::Agent, incidence::Incidence, age::Union{Integer, Nothing}=nothing,
    year::Union{Integer, Nothing}=nothing
)

    if isnothing(age)
        age = min(agent.age, incidence.max_age)
    end
    if isnothing(year)
        year = agent.cal_year
    end

    # assume no asthma if age < 3
    if age < 3
        has_asthma = false
    elseif age == 3
        has_asthma = rand(Bernoulli(prevalence_equation(
            agent.sex, age, year, agent.family_hist, agent.num_antibiotic_use, incidence
        )))
    else
        has_asthma = rand(Bernoulli(incidence_equation(
            agent.sex, age, year, agent.family_hist, agent.num_antibiotic_use, incidence
        )))
    end
    return has_asthma
end


function agent_has_asthma_prevalence(
    agent::Agent, incidence::Incidence, age::Union{Integer, Nothing}=nothing,
    year::Union{Integer, Nothing}=nothing
)

    if isnothing(age)
        age = min(agent.age - 1, incidence.max_age)
    end
    if isnothing(year)
        year = agent.cal_year - 1
    end

    # assume no asthma if age < 3
    if age < 3
        has_asthma = false
    else
        has_asthma = rand(Bernoulli(prevalence_equation(
            agent.sex, age, year, agent.family_hist, agent.num_antibiotic_use, incidence
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
function compute_asthma_age(agent::Agent, incidence::Incidence, current_age::Integer,
    max_asthma_age::Integer=110
)
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
        tmp_year = min(max(agent.cal_year - current_age + asthma_age, min_year), max_year)
        while find_asthma_age && asthma_age < max_asthma_age
            has_asthma = agent_has_asthma_incidence(agent, incidence, asthma_age, tmp_year)
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


function crude_incidence(sex, age, cal_year::Integer, parameters::AbstractDict)
    poly_age = poly_age_calculator(age)
    return exp(
        parameters[:β0] +
        parameters[:βsex] * sex +
        parameters[:βyear] * cal_year +
        parameters[:βsexyear] * sex * cal_year +
        sum(parameters[:βage] .* poly_age) +
        sum(parameters[:βsexage] .* sex .* poly_age)
    )
end


function incidence_equation(
    sex, age::Integer, cal_year::Integer, fam_hist::Bool, dose::Integer, incidence::Incidence
)
    parameters = incidence.parameters
    correction_year = min(cal_year, incidence.max_year + 1)
    cal_year = min(cal_year, incidence.max_year)
    p0 = crude_incidence(sex, age, cal_year, parameters)
    return inverse_logit(
        logit(p0) +
        fam_hist * log_OR_family_history(age, parameters[:βfam_hist]) +
        log_OR_abx_exposure(age, dose, parameters[:βabx_exp]) +
        incidence.incidence_correction_table[(correction_year, sex, min(age, 63))].correction[1]
    )
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


function prevalence_equation(
    sex, age::Integer, cal_year::Integer, family_hist::Bool, dose::Integer, incidence::Incidence
)
    parameters = incidence.parameters_prev
    correction_year = min(cal_year, incidence.max_year + 1)
    cal_year = min(cal_year, incidence.max_year)
    p0 = crude_prevalence(sex, age, cal_year, parameters)
    return inverse_logit(
        logit(p0) +
        family_hist * log_OR_family_history(age, parameters[:βfam_hist]) +
        log_OR_abx_exposure(age, dose, parameters[:βabx_exp]) +
        incidence.prevalence_correction_table[(correction_year, sex, min(age, 63))].correction[1]
    )
end


function crude_prevalence(sex, age::Integer, cal_year::Integer, parameters::AbstractDict)
    poly_year = poly_year_calculator(cal_year)
    poly_age = poly_age_calculator(age)
    poly_yearage = vec(poly_year .* poly_age')
    return exp(
        parameters[:β0] +
        parameters[:βsex] * sex +
        sum(parameters[:βyear] .* poly_year) +
        sum(parameters[:βage] .* poly_age) +
        sum(parameters[:βsexyear] .* sex .* poly_year) +
        sum(parameters[:βsexage] .* sex .* poly_age) +
        sum(parameters[:βyearage] .* poly_yearage) +
        sum(parameters[:βsexyearage] .* sex .* poly_yearage)
    )
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

function poly_year_calculator(
    year::Integer,
    alpha::Array{Float64, 1}=[2009.5, 2009.5],
    nd::Array{Float64, 1}=[1.0, 520.0, 17290.0, 456456.0]
)
    fs = zeros(3)
    fs[1] =  1 / sqrt(nd[2])
    fs[2] = (year-alpha[1]) / sqrt(nd[3])
    fs[3] = ((year-alpha[2]) * sqrt(nd[3]) * fs[2] - nd[3] / sqrt(nd[2]) * fs[1]) / sqrt(nd[4])
    popfirst!(fs)
    fs
end
