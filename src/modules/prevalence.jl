"""
    Prevalence

A struct containing information about asthma prevalence.

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
    `βyear::Vector{Float64}`: an array of 2 parameters to be multiplied by functions of year,
        i.e. βyear1 * g1(year) + βyear2 * g2(year)
        See `poly_year_calculator`.
    `βsexage::Vector{Float64}`: an array of 5 parameters to be multiplied by the sex and
        functions of age,
        i.e. βsexage1 * f1(age) * sex + βsexage2 * f2(age) * sex + βsexage3 * f3(age) * sex +
             βsexage4 * f4(age) * sex + βsexage5 * f5(age) * sex
        See `poly_age_calculator`.
    `βsexyear::Vector{Float64}`: an array of 2 parameters to be multiplied by sex and
        functions of year,
        i.e. βyear1 * g1(year) * sex + βyear2 * g2(year) * sex
    `βyearage::Vector{Float64}`: an array of 10 parameters to be multiplied by
        functions of age and year, i.e.
        βyearage1 * f1(age) * g1(year) + βyearage2 * f1(age) * g2(year) +
        βyearage3 * f2(age) * g1(year) + βyearage4 * f2(age) * g2(year) +
        βyearage5 * f3(age) * g1(year) + βyearage6 * f3(age) * g2(year) +
        βyearage7 * f4(age) * g1(year) + βyearage8 * f4(age) * g2(year) +
        βyearage9 * f5(age) * g1(year) + βyearage10 * f5(age) * g2(year)
    `βyearagesex::Vector{Float64}`: an array of 10 parameters to be multiplied by sex and
        functions of age and year, i.e.
        βyearagesex1 * f1(age) * g1(year) * sex + βyearagesex2 * f1(age) * g2(year) * sex +
        βyearagesex3 * f2(age) * g1(year) * sex + βyearagesex4 * f2(age) * g2(year) * sex +
        βyearagesex5 * f3(age) * g1(year) * sex + βyearagesex6 * f3(age) * g2(year) * sex +
        βyearagesex7 * f4(age) * g1(year) * sex + βyearagesex8 * f4(age) * g2(year) * sex +
        βyearagesex9 * f5(age) * g1(year) * sex + βyearagesex10 * f5(age) * g2(year) * sex
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
struct Prevalence <: PrevalenceModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    min_year::Integer
    max_year::Integer
    max_age::Integer
    correction_table::GroupedDataFrame{DataFrame}
    function Prevalence(config::AbstractDict)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:βage] = Array{Float64, 1}(parameters[:βage])
        parameters[:βyear] = Array{Float64, 1}(parameters[:βyear])
        parameters[:βsexage] = Array{Float64, 1}(parameters[:βsexage])
        parameters[:βsexyear] = Array{Float64, 1}(parameters[:βsexyear])
        parameters[:βyearage] = Array{Float64, 1}(parameters[:βyearage])
        parameters[:βsexyearage] = Array{Float64, 1}(parameters[:βsexyearage])
        parameters[:βfam_hist] = Array{Float64, 1}(parameters[:βfam_hist])
        parameters[:βabx_exp] = Array{Float64, 1}(parameters[:βabx_exp])
        max_age = config["max_age"]
        correction_table = load_prevalence_correction_table()
        min_year = collect(keys(correction_table)[1])[1]+1
        max_year = collect(keys(correction_table)[length(correction_table)])[1]
        new(
            hyperparameters, parameters, min_year, max_year, max_age,
            correction_table
        )
    end
    function Prevalence(hyperparameters::AbstractDict, parameters::AbstractDict,
        min_year::Integer, max_year::Integer, max_age::Integer,
        correction_table::GroupedDataFrame{DataFrame}
    )
        new(
            hyperparameters, parameters, min_year, max_year, max_age, correction_table
        )
    end
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
    agent_has_asthma(agent, prevalence)

Determine whether the agent obtains a new asthma diagnosis based on age and sex.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `prevalence::Prevalence`: Prevalence module, see [`Prevalence`](@ref).
"""
function agent_has_asthma(agent::Agent, prevalence::Prevalence)

    age = min(agent.age - 1, prevalence.max_age)
    year = agent.year - 1

    # assume no asthma if age < 3
    if age < 3
        has_asthma = false
    else
        has_asthma = rand(Bernoulli(prevalence_equation(
            agent.sex, age, year, agent.has_family_hist, agent.num_antibiotic_use, prevalence
        )))
    end
    return has_asthma
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
    sex, age::Integer, year::Integer, has_family_hist::Bool, dose::Integer,
    prevalence::Prevalence
)
    parameters = prevalence.parameters
    correction_year = min(year, prevalence.max_year + 1)
    year = min(year, prevalence.max_year)
    p0 = crude_prevalence(sex, age, year, parameters)
    p = logistic(
        logit(p0) +
        has_family_hist * log_OR_family_history(age, parameters[:βfam_hist]) +
        log_OR_abx_exposure(age, dose, parameters[:βabx_exp]) +
        prevalence.correction_table[(correction_year, sex, min(age, 63))].correction[1]
    )
    return p
end


function crude_prevalence(sex, age::Integer, year::Integer, parameters::AbstractDict)
    poly_year = poly_year_calculator(year)
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
