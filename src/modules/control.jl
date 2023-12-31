struct Control <: ControlModule
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
end


function process_control(agent::Agent, ctl::Control, initial::Bool=false)
    if initial
        age_scaled = (agent.age - 1) / 100
    else
        age_scaled = agent.age / 100
    end

    return control_prediction(
        (ctl.parameters[:β0] +
        age_scaled * ctl.parameters[:βage] +
        agent.sex*ctl.parameters[:βsex] +
        age_scaled * agent.sex * ctl.parameters[:βsexage] +
        age_scaled^2 * agent.sex * ctl.parameters[:βsexage2] +
        age_scaled^2 * ctl.parameters[:βage2]),
        ctl.parameters[:θ]
    )
end

# pred function
function control_prediction(eta::Float64, theta::Union{Float64,Vector{Float64}};
    inv_link::Function=StatsFuns.logistic)::Union{Float64,Vector{Float64}}
    theta = [-1e5;theta;1e5]
    [inv_link(theta[j+1] - eta) - inv_link(theta[j] - eta) for j in 1:(length(theta)-1)]
end


"""
    random_parameter_initialization!(control)

Assign the parameter β0 a random value from a normal distribution with a mean μ = β0_μ and a
standard deviation σ = β0_σ.

# Arguments

- `control:Control`: Control module, see [`Control`](@ref).
"""
function random_parameter_initialization!(control::Control)
    control.parameters[:β0] = rand(Normal(
        control.hyperparameters[:β0_μ],
        control.hyperparameters[:β0_σ]
    ))
end
