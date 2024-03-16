
struct AntibioticExposure <: AntibioticExposureModule
    hyperparameters::AbstractDict
    parameters::AbstractDict
    mid_trends::GroupedDataFrame{DataFrame}
    AbxOR
    function AntibioticExposure(config::AbstractDict)
        hyperparameters = string_to_symbols_dict(config["hyperparameters"])
        parameters = string_to_symbols_dict(config["parameters"])
        mid_trends = load_abx_mid_trends()
        new(hyperparameters, parameters, mid_trends, nothing)
    end
    function AntibioticExposure(hyperparameters::AbstractDict,
        parameters::AbstractDict, mid_trends::GroupedDataFrame, AbxOR)
        new(hyperparameters, parameters, mid_trends, AbxOR)
    end
end


function load_abx_mid_trends()
    abx_mid_trends = groupby(CSV.read(
        joinpath(PROCESSED_DATA_PATH, "midtrends.csv"),
        DataFrame
    ), [:year, :sex])
    return abx_mid_trends
end


function process_antibiotic_exposure(antibiotic_exposure::AntibioticExposure, sex::Bool, cal_year::Integer)
    if !isnothing(antibiotic_exposure.parameters[:fixyear])
        if isa(antibiotic_exposure.parameters[:fixyear], Number)
            return rand(
                NegativeBinomial(
                    antibiotic_exposure.parameters[:θ],
                    antibiotic_exposure_prob(
                        sex,
                        antibiotic_exposure.parameters[:fixyear],
                        antibiotic_exposure.parameters
            )))
        else
            tmp_mu = max(
                antibiotic_exposure.parameters[:midtrends][(cal_year, Int(sex))].rate[1],
                antibiotic_exposure.parameters[:βfloor]
            )
            return rand(NegativeBinomial(
                antibiotic_exposure.parameters[:θ],
                antibiotic_exposure.parameters[:θ] / (antibiotic_exposure.parameters[:θ] + tmp_mu)
            ))
        end
    else
        return rand(NegativeBinomial(
            antibiotic_exposure.parameters[:θ],
            antibiotic_exposure_prob(sex, cal_year, antibiotic_exposure.parameters)
        ))
    end

end


function process_antibiotic_exposure_initial(antibiotic_exposure::AntibioticExposure,
    sex::Bool, cal_born::Integer)

    if cal_born < 2001
        return rand(NegativeBinomial(
            antibiotic_exposure.parameters[:θ],
            antibiotic_exposure_prob(sex, 2000, antibiotic_exposure.parameters)
        ))
    else
        return process_antibiotic_exposure(antibiotic_exposure, sex, cal_born)
    end
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
    mu = exp(
        parameters[:β0] +
        parameters[:βsex] * sex +
        parameters[:βcal_year] * cal_year +
        parameters[:β2005] * (cal_year > 2005) +
        parameters[:β2005_cal_year] * (cal_year > 2005)*cal_year
    )
    mu = max(mu, parameters[:βfloor] / 1000)
    return parameters[:θ] / (parameters[:θ] + mu)
end

function random_parameter_initialization!(anti::AntibioticExposure)
    # anti.parameters[:β0] = rand(Normal(anti.hyperparameters[:β0_μ], anti.hyperparameters[:β0_σ]))
    nothing
end
