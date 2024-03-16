struct Exacerbation <: ExacerbationModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    initial_rate::Float64
    function Exacerbation(config::AbstractDict, province::String)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        initial_rate = config["initial_rate"]
        exacerbation_calibration = load_exacerbation_calibration()
        parameters[:calibration] = groupby(
            select(
                filter(
                    [:province] => (x) -> x == province,
                    exacerbation_calibration
                ),
                Not([:province])
            ),
            [:year,:sex]
        )
        parameters[:min_year] = collect(keys(parameters[:calibration])[1])[1]+1
        new(hyperparameters, parameters, initial_rate)
    end
    function Exacerbation(
        hyperparameters::AbstractDict, parameters::AbstractDict,
        initial_rate::Float64
    )
        new(hyperparameters, parameters, initial_rate)
    end
end

struct ExacerbationHist <: ExacerbationHistModule
    num_current_year::Integer
    num_prev_year::Integer
end

function load_exacerbation_calibration()
    exacerbation_calibration  = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_calibrated_exac.csv"),
        DataFrame
    )
    return exacerbation_calibration
end


"""
    compute_num_exacerbations(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(agent::Agent, exacerbation::Exacerbation)
    return compute_num_exacerbations(agent.age, agent.sex, agent.cal_year, agent.control_levels,
        exacerbation::Exacerbation)
end


"""
    compute_num_exacerbations(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(age::Integer, sex::Bool, cal_year::Integer,
    control_levels::AbstractDict, exacerbation::Exacerbation)
    params = exacerbation.parameters
    tmp_year = max(params[:min_year], cal_year)
    tmp_age = min(age, 90)
    μ = (
        params[:β0] +
        params[:β0_calibration] +
        age * params[:βage] +
        sex * params[:βsex] +
        control_levels[:uncontrolled] * params[:βcontrol_UC] +
        control_levels[:partially_controlled] * params[:βcontrol_PC] +
        control_levels[:fully_controlled] * params[:βcontrol_C] +
        log(params[:calibration][(tmp_year, Int(sex))][tmp_age - 2, "calibrator_multiplier"])
    )
    return exacerbation_prediction(μ)
end


"""
    compute_num_exacerbations_initial(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations_initial(agent::Agent, exacerbation::Exacerbation)
    params = exacerbation.parameters
    tmp_year = max(params[:min_year], agent.cal_year - 1)
    tmp_age = min(agent.age - 1, 90)
    if tmp_age < 3
        return 0
    else
        μ = (
            params[:β0] +
            tmp_age * params[:βage] +
            agent.sex * params[:βsex] +
            agent.control_levels[:uncontrolled] * params[:βcontrol_UC] +
            agent.control_levels[:partially_controlled] * params[:βcontrol_PC] +
            agent.control_levels[:fully_controlled] * params[:βcontrol_C] +
            log(params[:calibration][(tmp_year, Int(agent.sex))][tmp_age - 2, "calibrator_multiplier"])
        )
        return exacerbation_prediction(μ)
    end
end

function random_parameter_initialization!(exac::Exacerbation)
    exac.parameters[:β0] = rand(Normal(exac.hyperparameters[:β0_μ], exac.hyperparameters[:β0_σ]))
end

function exacerbation_prediction(μ::Float64; inv_link::Function=exp)
    rand(Poisson(inv_link(μ)))
end
