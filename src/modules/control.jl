struct Control <: ControlModule
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
end


function process_control(control::Control, sex::Bool, age::Integer, initial::Bool=false)
    if initial
        age_scaled = (age - 1) / 100
    else
        age_scaled = age / 100
    end

    return control_prediction(
        (
            control.parameters[:β0] +
            age_scaled * control.parameters[:βage] +
            sex * control.parameters[:βsex] +
            age_scaled * sex * control.parameters[:βsexage] +
            age_scaled^2 * sex * control.parameters[:βsexage2] +
            age_scaled^2 * control.parameters[:βage2]
        ),
        control.parameters[:θ]
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
