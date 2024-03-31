"""
    AntibioticExposure

A struct containing information about antibiotic use.

# Fields
- `parameters::AbstractDict`: A dictionary with the following keys:
    `β0`: float, the constant parameter when computing μ.
    `βcal_year`: float, the parameter to be multiplied by the agent's birth year for computing μ.
    `β2005`: float, an added constant parameter if the agent's birth year > 2005 for computing μ.
    `βsex`: float, the parameter to be multiplied by the agent's sex when computing μ.
    `θ`: int, the number of successes (the r parameter) in the negative binomial distribution.
    `β2005_cal_year`: float, if the agent's birth year is > 2005, β2005_cal_year will be
        multiplied by the birth year when computing μ.
    `fixyear`: integer or nothing. If present, replaces the `year` parameter when computing
        the probability for the negative binomial distribution.
    `βfloor`: float, the minimum value of μ.
- `mid_trends::GroupedDataFrame{DataFrame}`: a set of data frames grouped by year and sex.
    Each entry is a DataFrame with a single row with the following columns:
    `year`: integer
    `sex`: 0 = female, 1 = male
    `rate`: float, TODO.
"""
struct AntibioticExposure <: AntibioticExposureModule
    parameters::AbstractDict
    mid_trends::GroupedDataFrame{DataFrame}
    function AntibioticExposure(config::AbstractDict)
        parameters = string_to_symbols_dict(config["parameters"])
        mid_trends = load_abx_mid_trends()
        new(parameters, mid_trends)
    end
    function AntibioticExposure(parameters::AbstractDict, mid_trends::GroupedDataFrame)
        new(parameters, mid_trends)
    end
end


"""
    load_abx_mid_trends()

Load the antibiotic mid trends table.

# Returns
- `GroupedDataFrame{DataFrame}`: a set of data frames grouped by year and sex. Each entry is a
    DataFrame with a single row with the following columns:
    `year`: integer
    `sex`: 0 = female, 1 = male
    `rate`: float, TODO.
"""
function load_abx_mid_trends()
    abx_mid_trends = groupby(CSV.read(
        joinpath(PROCESSED_DATA_PATH, "midtrends.csv"),
        DataFrame
    ), [:year, :sex])
    return abx_mid_trends
end


"""
    compute_num_antibiotic_use(antibiotic_exposure, sex, birth_year)

Compute the number of antibiotics used during the first year of life.

# Arguments
- `antibiotic_exposure::AntibioticExposure`: contains information about antibiotic exposure,
    see [`AntibioticExposure`](@ref).
- `sex::Bool`: Sex of agent, true = male, false = female.
- `birth_year::Integer`: The year the agent was born.
"""
function compute_num_antibiotic_use(
    ; antibiotic_exposure::AntibioticExposure, sex::Bool, birth_year::Integer
)
    if birth_year < 2001
        p = antibiotic_exposure_prob(
            sex=sex, cal_year=2000, parameters=antibiotic_exposure.parameters
        )
    elseif !isnothing(antibiotic_exposure.parameters[:fixyear])
        if isa(antibiotic_exposure.parameters[:fixyear], Number)
            p = antibiotic_exposure_prob(
                sex=sex,
                cal_year=antibiotic_exposure.parameters[:fixyear],
                parameters=antibiotic_exposure.parameters
            )
        else
            μ = max(
                antibiotic_exposure.parameters[:midtrends][(birth_year, Int(sex))].rate[1],
                antibiotic_exposure.parameters[:βfloor]
            )
            p = antibiotic_exposure.parameters[:θ] / (antibiotic_exposure.parameters[:θ] + μ)
        end
    else
        p = antibiotic_exposure_prob(
            sex=sex,
            cal_year=birth_year,
            parameters=antibiotic_exposure.parameters
        )
    end
    r = antibiotic_exposure.parameters[:θ]
    return rand(NegativeBinomial(r, p))
end


"""
    antibiotic_exposure_prob(sex, cal_year, parameters)

Returns the probability of antibiotic exposure for a given year and sex.

# Arguments
- `sex::Bool`: Sex of agent, true = male, false = female.
- `cal_year::Integer`: The calendar year.
- `parameters::AbstractDict`: A dictionary with the following keys:
    `β0`: float, the constant parameter when computing μ.
    `βcal_year`: float, the parameter to be multiplied by the agent's birth year for computing μ.
    `β2005`: float, an added constant parameter if the agent's birth year > 2005 for computing μ.
    `βsex`: float, the parameter to be multiplied by the agent's sex when computing μ.
    `θ`: int, the number of successes (the r parameter) in the negative binomial distribution.
    `β2005_cal_year`: float, if the agent's birth year is > 2005, β2005_cal_year will be
        multiplied by the birth year when computing μ.
    `fixyear`: integer or nothing. If present, replaces the `year` parameter when computing
        the probability for the negative binomial distribution.
    `βfloor`: float, the minimum value of μ.
"""
function antibiotic_exposure_prob(; sex::Bool, cal_year::Integer, parameters::AbstractDict)
    μ = exp(
        parameters[:β0] +
        parameters[:βsex] * sex +
        parameters[:βcal_year] * cal_year +
        parameters[:β2005] * (cal_year > 2005) +
        parameters[:β2005_cal_year] * (cal_year > 2005) * cal_year
    )
    μ = max(μ, parameters[:βfloor] / 1000)
    return parameters[:θ] / (parameters[:θ] + μ)
end
