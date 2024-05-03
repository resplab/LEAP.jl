"""
    Exacerbation

A struct containing information about asthma exacerbations.

# Fields
- `hyperparameters::AbstractDict`: A dictionary containing the hyperparameters used
    to compute `β0` from a normal distribution:
    `β0_μ`: Float64, the mean of the normal distribution.
    `β0_σ`: Float64, the standard deviation of the normal distribution.
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `β0`: Float64, a constant parameter, randomly selected from a normal distribution with mean
        `β0_μ` and standard deviation `β0_σ`. See `hyperparameters`.
    `β0_calibration`: Float64, the parameter for the calibration term.
    `βage`: Float64, the parameter for the age term.
    `βsex`: Float64, the parameter for the sex term.
    `βasthmaDx`: Float64, TODO.
    `βprev_exac1`: Float64, TODO.
    `βprev_exac2`: Float64, TODO.
    `βcontrol`: Float64, TODO.
    `βcontrol_C`: Float64, the parameter for the controlled asthma term.
    `βcontrol_PC`: Float64, the parameter for the partially-controlled asthma term.
    `βcontrol_UC`: Float64, the parameter for the uncontrolled asthma term.
    `min_year`: Integer, the minimum year for which exacerbation data exists + 1. Currently 2001.
- `calibration_table::GroupedDataFrame{DataFrame}`: A dataframe grouped by year and sex, with the
    following columns:
        `year`: integer year.
        `sex`: 1 = male, 0 = female.
        `age`: integer age.
        `calibrator_multiplier`: Float64, TODO.
    See `master_calibrated_exac.csv`.
- `inital_rate::Float64`: TODO.
"""
struct Exacerbation <: ExacerbationModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    calibration_table::GroupedDataFrame{DataFrame}
    initial_rate::Float64
    function Exacerbation(config::AbstractDict, province::String)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:β0] = assign_random_β0(hyperparameters[:β0_μ], hyperparameters[:β0_σ])
        initial_rate = config["initial_rate"]
        calibration_table = load_exacerbation_calibration(province)
        parameters[:min_year] = collect(keys(calibration_table)[1])[1] + 1
        new(hyperparameters, parameters, calibration_table, initial_rate)
    end
    function Exacerbation(
        hyperparameters::AbstractDict, parameters::AbstractDict,
        calibration_table::GroupedDataFrame{DataFrame}, initial_rate::Float64
    )
        new(hyperparameters, parameters, calibration_table, initial_rate)
    end
end

struct ExacerbationHist <: ExacerbationHistModule
    num_current_year::Integer
    num_prev_year::Integer
end


"""
    load_exacerbation_calibration(province)

Load the exacerbation calibration table.

# Arguments
- `province::String`: a string indicating the province abbreviation, e.g. "BC". For all of Canada,
    set province to "CA".

# Returns
- `GroupedDataFrame{DataFrame}`: A dataframe grouped by year and sex, with the
    following columns:
    `year`: integer year.
    `sex`: 1 = male, 0 = female.
    `age`: integer age.
    `calibrator_multiplier`: Float64, TODO.
"""
function load_exacerbation_calibration(province::String)
    exacerbation_calibration  = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_calibrated_exac.csv"),
        DataFrame
    )
    table = groupby(
        select(
            filter(
                [:province] => (x) -> x == province,
                exacerbation_calibration
            ),
            Not([:province])
        ),
        [:year, :sex]
    )
    return table
end


"""
    compute_num_exacerbations(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(agent::AgentModule, exacerbation::Exacerbation)
    return compute_num_exacerbations(agent.age, agent.sex, agent.cal_year, agent.control_levels,
        exacerbation::Exacerbation)
end


"""
    compute_num_exacerbations(age, sex, cal_year, control_levels, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(age::Integer, sex::Bool, cal_year::Integer,
    control_levels::AbstractDict, exacerbation::Exacerbation)
    params = exacerbation.parameters
    year = max(params[:min_year], cal_year)
    tmp_age = min(age, 90)
    μ = (
        params[:β0] +
        params[:β0_calibration] +
        age * params[:βage] +
        sex * params[:βsex] +
        control_levels[:uncontrolled] * params[:βcontrol_UC] +
        control_levels[:partially_controlled] * params[:βcontrol_PC] +
        control_levels[:fully_controlled] * params[:βcontrol_C] +
        log(exacerbation.calibration_table[(year, Int(sex))][tmp_age - 2, "calibrator_multiplier"])
    )
    λ = exp(μ)
    poisson_distribution = Poisson(λ)
    return rand(poisson_distribution)
end


"""
    compute_num_exacerbations_initial(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations_initial(agent::AgentModule, exacerbation::Exacerbation)
    params = exacerbation.parameters
    year = max(params[:min_year], agent.cal_year - 1)
    age = min(agent.age - 1, 90)
    if age < 3
        return 0
    else
        μ = (
            params[:β0] +
            age * params[:βage] +
            agent.sex * params[:βsex] +
            agent.control_levels[:uncontrolled] * params[:βcontrol_UC] +
            agent.control_levels[:partially_controlled] * params[:βcontrol_PC] +
            agent.control_levels[:fully_controlled] * params[:βcontrol_C] +
            log(exacerbation.calibration_table[(year, Int(agent.sex))][age - 2, "calibrator_multiplier"])
        )
        λ = exp(μ)
        poisson_distribution = Poisson(λ)
        return rand(poisson_distribution)
    end
end


"""
    assign_random_β0(β0_μ, β0_σ)

Assign the parameter β0 a random value from a normal distribution with a mean μ = β0_μ and a
standard deviation σ = β0_σ.

# Arguments

- `β0_μ::Float64`: The mean for the normal distribution.
- `β0_σ::Float64`: The standard deviation for the normal distribution.
"""
function assign_random_β0(β0_μ::Float64, β0_σ::Float64)
    β0 = rand(Normal(β0_μ, β0_σ))
    return β0
end

function assign_random_β0!(exacerbation::Exacerbation)
    β0 = assign_random_β0(
        exacerbation.hyperparameters[:β0_μ],
        exacerbation.hyperparameters[:β0_σ]
    )
    exacerbation.parameters[:β0] = β0
end
