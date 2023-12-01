struct Exacerbation <: Exacerbation_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    initial_rate::Float64
end

"""
    compute_num_exacerbations(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(agent::Agent, exacerbation::Exacerbation)
    return compute_num_exacerbations(agent.age, agent.sex, agent.cal_year, agent.control,
        exacerbation::Exacerbation)
end


"""
    compute_num_exacerbations(agent, exacerbation)

Return the number of asthma exacerbations in a given year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exacerbation::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function compute_num_exacerbations(age::Integer, sex::Bool, cal_year::Integer, control::Control
    exacerbation::Exacerbation)
    params = exacerbation.parameters
    tmp_year = max(params[:min_year], cal_year)
    tmp_age = min(age, 90)
    μ = (
        params[:β0] +
        params[:β0_calibration]
        age * params[:βage] +
        sex * params[:βsex] +
        control[1] * params[:βcontrol_C] +
        control[2] * params[:βcontrol_PC] +
        control[3] * params[:βcontrol_UC] +
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
            agent.control[1] * params[:βcontrol_C] +
            agent.control[2] * params[:βcontrol_PC] +
            agent.control[3] * params[:βcontrol_UC] +
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
