# include("abstractModule.jl")
# include("../utils")

"""
    AntibioticExposure

A struct for antibiotic exposure.

# Fields
- `hyperparameters::Union{AbstractDict, Nothing}`: A dictionary containing two keys,
    `β0_μ` and `β0_σ`.
- `parameters::Union{AbstractDict,Nothing}`: A dictionary containing the following keys:
    `θ`, `β0`, `βage`, `βsex`, `βcal_year`, `β2005`, `β2005_cal_year`, `fixyear`, `fix2000`,
    `βfloor`, `midtrends`.
- `AbxOR: TODO.
"""
struct AntibioticExposure <: AntibioticExposure_Module
    hyperparameters::AbstractDict
    parameters::AbstractDict
    AbxOR
end


"""
    process(agent, antibio)

Return a random value from a negative binomial distribution.

# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `antibio::AntibioticExposure`: AntibioticExposure module, see [`AntibioticExposure`](@ref).

"""
function process(agent::Agent, antibio::AntibioticExposure)
    if !isnothing(antibio.parameters[:fixyear])
        if isa(antibio.parameters[:fixyear], Number)
            r = antibio.parameters[:θ]
            p = antibiotic_exposure_prob(
                agent.sex, antibio.parameters[:fixyear], antibio.parameters
            )
        else
            tmp_mu = max(
                antibio.parameters[:midtrends][(agent.cal_year, Int(agent.sex))].rate[1],
                antibio.parameters[:βfloor]
            )
            r = anti.parameters[:θ]
            p = anti.parameters[:θ] / (antibio.parameters[:θ] + tmp_mu)
        end
    else
        r = anti.parameters[:θ]
        p = antibiotic_exposure_prob(agent.sex, agent.cal_year, antibio.parameters)
    end
    return rand(NegativeBinomial(r, p))
end


"""
    process_initial(agent, antibio)

Return a random value from a negative binomial distribution.

# Arguments
- `agent::Agent`: A person in the model, see [`Agent`](@ref).
- `antibio::AntibioticExposure`: AntibioticExposure module, see [`AntibioticExposure`](@ref).

"""
function process_initial(agent::Agent, antibio::AntibioticExposure, cal_born::Integer)
    if cal_born < 2001
        r = antibio.parameters[:θ]
        p = antibiotic_exposure_prob(
            agent.sex, 2000, antibio.parameters
        )
    else
        if !isnothing(antibio.parameters[:fixyear])
            if isa(antibio.parameters[:fixyear], Number)
                r = antibio.parameters[:θ]
                p = antibiotic_exposure_prob(
                    agent.sex, antibio.parameters[:fixyear], antibio.parameters
                )
            else
                tmp_mu = max(
                    antibio.parameters[:midtrends][(cal_born, Int(agent.sex))].rate[1],
                    antibio.parameters[:βfloor]
                )
                r = anti.parameters[:θ]
                p = anti.parameters[:θ] / (antibio.parameters[:θ] + tmp_mu)
            end
        else
            r = anti.parameters[:θ]
            p = antibiotic_exposure_prob(agent.sex, cal_born, antibio.parameters)
        end
    end
    return rand(NegativeBinomial(r, p))
end


"""
    antibiotic_exposure_prob(sex, cal_year, parameters)

Returns the probability of antibiotic exposure for a given year and sex.

# Arguments
- `sex::Bool`: Sex of agent, true = male, false = female.
- `cal_year::Integer`: The calendar year.
- `parameters::AbstractDict`: A dictionary containing the following keys:
    `θ`, `β0`, `βage`, `βsex`, `βcal_year`.

"""
function antibiotic_exposure_prob(sex::Bool, cal_year::Integer, parameters::AbstractDict)
    μ = (
        exp(parameters[:β0] +
        parameters[:βsex] * sex +
        parameters[:βcal_year] * cal_year +
        parameters[:β2005] * (cal_year > 2005) +
        parameters[:β2005_cal_year] * (cal_year > 2005) * cal_year
    )
    μ = max(μ, parameters[:βfloor] / 1000)
    return parameters[:θ] / (parameters[:θ] + μ)
end

function random_parameter_initialization!(antibio::AntibioticExposure)
    # anti.parameters[:β0] = rand(Normal(anti.hyperparameters[:β0_μ], anti.hyperparameters[:β0_σ]))
    nothing
end
