
"""
    Exacerbation

A struct containing information about asthma exacerbations.

# Fields
- `hyperparameters::Union{AbstractDict,Nothing}`: A dictionary containing two keys,
    `β0_μ` and `β0_σ`.
- `parameters::Union{AbstractDict,Nothing}`: A dictionary containing the following keys:
    `β0`, `βage`, `βsex`, `βasthmaDx`, `βprev_exac1,` `βprev_exac2`, `βcontrol_C`, `βcontrol_PC`,
    and `βcontrol_UC`.
- `initial_rate::Float64`: Initial asthma exacerbation rate per year.
"""
struct Exacerbation <: Exacerbation_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    initial_rate::Float64
end



"""
    process(agent, exacerbation, inv_link)

Given an agent and asthma exacerbation details, assign the exacerbation rate per year.


# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `exac::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).
"""
function process(agent::Agent, exac::Exacerbation)
    parameters = exac.parameters
    tmp_cal_year = max(agent.cal_year, parameters[:min_year])
    tmp_age = min(agent.age, 90)
    μ = (
        parameters[:β0] +
        agent.age * parameters[:βage] +
        agent.sex * parameters[:βsex] +
        parameters[:β0_calibration] +
        agent.control[1] * parameters[:βcontrol_C] +
        agent.control[2] * parameters[:βcontrol_PC] +
        agent.control[3] * parameters[:βcontrol_UC] +
        log(parameters[:calibration][(tmp_cal_year, Int(ag.sex))][tmp_age-2,"calibrator_multiplier"]))
    )
    return exacerbation_prediction(μ)
end

"""
    process_initial(agent, exacerbation)

Given an agent and asthma exacerbation details, assign the exacerbation rate per year.

# Arguments
- `agent::Agent`: A person in the model, see  [`Agent`](@ref).
- `exac::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).

"""
function process_initial(agent::Agent, exac::Exacerbation)
    parameters = exac.parameters
    tmp_year = max(parameters[:min_year], agent.cal_year-1)
    tmp_age = min(agent.age-1,90)
    if tmp_age <3
        return 0
    else
        μ = (
            parameters[:β0] +
            agent.age * parameters[:βage] +
            agent.sex * parameters[:βsex] +
            parameters[:β0_calibration] +
            agent.control[1] * parameters[:βcontrol_C] +
            agent.control[2] * parameters[:βcontrol_PC] +
            agent.control[3] * parameters[:βcontrol_UC] +
            log(parameters[:calibration][(tmp_year, Int(agent.sex))][tmp_age-2,"calibrator_multiplier"]))
        )
        return exacerbation_prediction(μ)
    end
end


"""
    random_parameter_initialization!(exac)

Assign the parameter β0 a random value from a normal distribution with a mean μ = β0_μ and a
standard deviation σ = β0_σ.

# Arguments
- `exac::Exacerbation`: An asthma exacerbation, see [`Exacerbation`](@ref).

"""
function random_parameter_initialization!(exac::Exacerbation)
    exac.parameters[:β0] = rand(
        Normal(μ=exac.hyperparameters[:β0_μ], σ=exac.hyperparameters[:β0_σ])
    )
end


"""
    exacerbation_prediction(μ, inv_link)

Assign the exacerbation rate per year.


# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `inv_link::Function=exp`: A function to apply, default is exponential.
"""
function exacerbation_prediction(μ::Float64; inv_link::Function=exp)
    return rand(Poisson(inv_link(μ)))
end
