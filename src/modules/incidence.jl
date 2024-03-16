"""
    Incidence

A struct containing information about asthma incidence.

# Fields
- `hyperparameters::AbstractDict`: A dictionary containing two keys,
    `β0_μ` and `β0_σ`.
- `parameters::Union{AbstractDict,Nothing}`: A dictionary containing the following keys:
    `βage`, `βage2`,`βage3`, `βage4`, `βage5`, `βageM`, `βage2M`, `βage3M`, `βage4M`, `βage5M`,
    `βsex`, `matrix_Abx_risk`, `βcal_year`, `βcal_yearM`, `βCABE`, `β0`, `β0_correction`,
    and `β0_overall_correction`.
- `incidence_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
    by year and province. Each data frame contains the following columns:
        `year`: integer year.
        `age`: integer age.
        `F`: Float64.
        `M`: Float64.
        `province`: A string indicating the province abbreviation, e.g. "BC".
    See `master_incidence_rate`.
- `prevalence_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
    by year and province. Each data frame contains the following columns:
        `year`: integer year.
        `age`: integer age.
        `F`: Float64.
        `M`: Float64.
        `province`: A string indicating the province abbreviation, e.g. "BC".
    See `master_prevalence_rate`.
- `calibration_table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped
    by year, sex, fam_history, abx_exposure. Each data frame contains the following columns:
        `fam_history`: integer of 0 or 1, 0 = no family history, 1 = family history.
        `abx_exposure`: integer indicating the number of antibiotic prescriptions during
            the first year of life.
        `year`: integer year.
        `sex`: sex of person, 1 = male, 0 = female.
        `age`: integer age.
        `prob`: Float64.
        `OR`: Float64.
        `calibrated_prev`: Float64, the calibrated asthma prevalence.
        `prev`: Float64.
        `calibrated_inc`: Float64, the calibrated asthma incidence.
        `province`: A string indicating the province abbreviation, e.g. "BC".
    See `M3_calibrated_asthma_prev_inc`.
"""
struct Incidence <: IncidenceModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    incidence_table::GroupedDataFrame{DataFrame}
    prevalence_table::GroupedDataFrame{DataFrame}
    calibration_table::GroupedDataFrame{DataFrame}
    min_year::Integer
    max_year::Integer
    initial_distribution
    function Incidence(config::AbstractDict, starting_year::Integer, province::String)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        incidence_table = load_incidence_table(starting_year, province)
        prevalence_table = load_prevalence_table(starting_year, province)
        calibration_table = load_calibration_table(province)
        min_year = collect(keys(calibration_table)[1])[1]+1
        max_year = collect(keys(calibration_table)[length(calibration_table)])[1]
        new(hyperparameters, parameters, incidence_table, prevalence_table, calibration_table,
            min_year, max_year, nothing)
    end
    function Incidence(hyperparameters::AbstractDict, parameters::AbstractDict,
        incidence_table::GroupedDataFrame{DataFrame}, prevalence_table::GroupedDataFrame{DataFrame},
        calibration_table::GroupedDataFrame{DataFrame}, min_year::Integer,
        max_year::Integer, initial_distribution)
        new(hyperparameters, parameters, incidence_table, prevalence_table, calibration_table,
            min_year, max_year, initial_distribution)
    end
end


function load_incidence_table(starting_year::Integer, province::String)
    master_incidence_rate = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_inc_interpolated.csv"),
        DataFrame
    )
    incidence_table = groupby(
        filter([:year,:province] => (x,y) -> x >= min(
            starting_year, master_incidence_rate.year[nrow(master_incidence_rate)]
            ) && y==province, master_incidence_rate
            ),
        :year
    )
    return incidence_table
end


function load_prevalence_table(starting_year::Integer, province::String)
    master_prevalence_rate = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_prev_interpolated.csv"),
        DataFrame
    )
    prevalence_table = groupby(
        filter([:year, :province] => (x,y) -> x >= min(
            starting_year-1, master_prevalence_rate.year[nrow(master_prevalence_rate)]
            ) && y==province, master_prevalence_rate),
        :year
    )
    return prevalence_table
end


function load_calibration_table(province::String)
    M3_calibrated_asthma_prev_inc = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_calibrated_asthma_prev_inc_M3.csv"),
        DataFrame
    )
    calibration_table = groupby(
        select(
            filter([:province] => (x) -> x == province, M3_calibrated_asthma_prev_inc),
            Not([:province])
        ),
        [:year, :sex, :fam_history, :abx_exposure]
    )
    return calibration_table
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
        inc_or_prev_cal = "calibrated_prev"
    elseif inc_or_prev == "incidence"
        inc_or_prev_cal = "calibrated_inc"
    else
        throw(ArgumentError(
            "inc_or_prev must be either 'prevalence' or 'incidence', received $inc_or_prev."
        ))
    end

    # assume no asthma if age < 3
    if tmp_age < 3
        has_asthma = false
    elseif tmp_age < 7
        has_asthma = rand(Bernoulli(incidence.calibration_table[(
            tmp_year, Int(agent.sex), Int(agent.family_hist), min(agent.num_antibiotic_use, 3)
            )][tmp_age - 2, inc_or_prev_cal]))
    else # no effect of Abx beyond 7 years of age
        has_asthma = rand(Bernoulli(incidence.calibration_table[(
            tmp_year, Int(agent.sex), Int(agent.family_hist), 0
            )][tmp_age - 2, inc_or_prev_cal]))
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
                incidence.calibration_table[(
                    tmp_year, tmp_sex, tmp_family_hist,
                    (asthma_age < 7 ? min(tmp_abx_num,3) : 0)
                )][asthma_age-2,"calibrated_inc"]))
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

function incidence_logit_prob(p::Real, sex::Bool, age::Int64, cal_year::Int64, ABX::Int64,
    parameters::AbstractDict, abxlogOR, calibration)
    return p
end


function random_parameter_initialization!(inc::Incidence)
    inc.parameters[:β0] = rand(Normal(inc.hyperparameters[:β0_μ], inc.hyperparameters[:β0_σ]))
end
