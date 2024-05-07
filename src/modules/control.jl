"""
    Control

A struct containing information about asthma control. This refers to how well the condition is
managed.
There are three levels of asthma control:
    fully-controlled = 1
    partially-controlled = 2
    uncontrolled = 3

# Fields
- `hyperparameters::AbstractDict`: A dictionary containing the hyperparameters used
    to compute `β0` from a normal distribution:
    `β0_μ`: Float64, the mean of the normal distribution.
    `β0_σ`: Float64, the standard deviation of the normal distribution.
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `β0`: Float64, a constant parameter. See `hyperparameters`.
    `βage`: Float64, the parameter for the age term.
    `βsex`: Float64, the parameter for the sex term.
    `βsexage`: Float64, the parameter for the sex * age term.
    `βsexage2`: Float64, the parameter for the sex * age^2 term.
    `βage2`: Float64, the parameter for the age^2 term.
    `βDx2`: Float64, unused?
    `βDx3`: Float64, unused?
    `θ`: An array of two numbers, which are used as the thresholds to compute the ordinal
        regression.
"""
@kwdef struct Control <: ControlModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    function Control(config::AbstractDict)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        parameters[:θ] = Array{Float64, 1}(parameters[:θ])
        parameters[:β0] = assign_random_β0(hyperparameters[:β0_μ], hyperparameters[:β0_σ])
        new(hyperparameters, parameters)
    end
    function Control(hyperparameters::AbstractDict, parameters::AbstractDict)
        new(hyperparameters, parameters)
    end
end


"""
    compute_control_levels(control, sex, age, initial)

Compute the probability that the control level = k for each value of k.
The probability is given by ordinal regression, where y = control level:
    P(y <= k) = σ(θ_k - η)
    P(y == k) = P(y <= k) - P(y < k + 1)
              = σ(θ_k - η) - σ(θ_(k+1) - η)

# Arguments
- `control:Control`: Control module, see [`Control`](@ref).
- `sex::Bool`: Sex of person, true = male, false = female.
- `age::Integer`: The age of the person (agent) in years.
- `initial::Bool`: if this is the initial computation.

# Returns
- `AbstractDict`: a dictionary with the probability of each control level.
    For example:
    {
        "fully_controlled": 0.2,
        "partially_controlled": 0.75,
        "uncontrolled": 0.05
    }
"""
function compute_control_levels(control::Control, sex::Bool, age::Integer,
    initial::Bool=false)::AbstractDict

    if initial
        age_scaled = (age - 1) / 100
    else
        age_scaled = age / 100
    end

    η = (
        control.parameters[:β0] +
        age_scaled * control.parameters[:βage] +
        sex * control.parameters[:βsex] +
        age_scaled * sex * control.parameters[:βsexage] +
        age_scaled^2 * sex * control.parameters[:βsexage2] +
        age_scaled^2 * control.parameters[:βage2]
    )
    control_levels_prob = compute_ordinal_regression(η, control.parameters[:θ])
    control_levels = Dict(
        :fully_controlled => control_levels_prob[1],
        :partially_controlled => control_levels_prob[2],
        :uncontrolled => control_levels_prob[3],
        :as_array => control_levels_prob
    )
    return control_levels
end


"""
    compute_ordinal_regression(η, θ, prob_function)

Compute the probability that y = k for each value of k.
The probability is given by ordinal regression:
    P(y <= k) = σ(θ_k - η)
    P(y == k) = P(y <= k) - P(y < k + 1)
              = σ(θ_k - η) - σ(θ_(k+1) - η)

# Arguments

- `η::Float64`: the weight for the regression.
- `θ::Union{Float64, Vector{Float64}}`: either a single value or an array of values for the
    threshold parameter.
- `prob_function::Function: A function to apply, default is the StatsFuns.logistic sigmoid function.

# Returns
- `Vector{Float64}`: a vector with the probability of each value of k.
    For example:

    k=1 | k=2  | k=2
    0.2 | 0.75 | 0.05
"""
function compute_ordinal_regression(η::Float64, θ::Union{Float64, Vector{Float64}};
    prob_function::Function=StatsFuns.logistic)::Vector{Float64}

    θ = [-1e5; θ; 1e5]
    return [
        prob_function(θ[k + 1] - η) -
        prob_function(θ[k] - η)
        for k in 1:(length(θ) - 1)
    ]
end


"""
    assign_random_β0(β0_μ, β0_σ)

Assign the parameter β0 a random value from a normal distribution with a mean μ = β0_μ and a
standard deviation σ = β0_σ.

# Arguments

- `β0_μ::Float64`: The mean for the normal distribution.
- `β0_σ::Float64`: The standard deviation for the normal distribution.
"""
function assign_random_β0(β0_μ::Float64, β0_σ::Float64)::Float64
    β0 = rand(Normal(β0_μ, β0_σ))
    return β0
end


function assign_random_β0!(control::Control)
    β0 = assign_random_β0(
        control.hyperparameters[:β0_μ],
        control.hyperparameters[:β0_σ]
    )
    control.parameters[:β0] = β0
end
