import Base.@kwdef


"""
    Agent

A person in the model.

# Fields
- `sex::Bool`: Sex of person, true = male, false = female.
- `age::Integer`: Age of person in years.
- `cal_year::Integer`: the calendar year of the current iteration, e.g. 2027.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if
    the simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.
- `alive::Bool`: Whether the person is alive, true = alive.
- `num_antibiotic_use::Int`: TODO.
- `has_asthma::Bool`: Whether the person has astham, true = has asthma.
- `asthma_age::Union{Nothing, Int}`: Age at which the person was diagnosed with asthma.
- `severity::Union{Nothing, Int}`: Asthma severity level: 1 = mild, 2 = severe, 3 = very severe.
- `control_levels::Union{Nothing, AbstractDict}`: Asthma control level: 1 = fully controlled,
    2 = partially controlled, 3 = uncontrolled.
- `exac_hist::Union{Nothing, ExacerbationHistModule}`: Total number of exacerbations.
- `exac_sev_hist::Union{Nothing, ExacerbationSeverityHistModule}`: Number of exacerbations by severity.
- `total_hosp::Integer`: Total number of very severe asthma exacerbations leading to hospitalization.
- `has_family_hist::Bool`: Is there a family history of asthma?
- `asthma_status::Bool`: TODO.
- `census_division::CensusDivisionModule`: the Canadian census division where the agent resides.
"""
@kwdef struct Agent  <: AgentModule
    sex::Bool
    age::Integer
    cal_year::Integer
    cal_year_index::Integer
    alive::Bool
    num_antibiotic_use::Integer
    has_asthma::Bool
    asthma_age::Union{Nothing, Integer}
    severity::Union{Nothing, Integer}
    control_levels::Union{Nothing, AbstractDict}
    exac_hist::Union{Nothing, ExacerbationHistModule}
    exac_sev_hist::Union{Nothing, ExacerbationSeverityHistModule}
    total_hosp::Integer
    has_family_hist::Bool
    asthma_status::Bool
    census_division::Union{Nothing, CensusDivisionModule}
    pollution::Union{Nothing, PollutionModule}
    function Agent(; sex::Bool, age::Integer, cal_year::Integer,
        cal_year_index::Integer, alive::Bool=true, has_asthma::Bool=false,
        asthma_age::Union{Integer, Nothing}=nothing, severity=nothing, control_levels=nothing,
        exac_hist::ExacerbationHistModule=ExacerbationHist(0, 0),
        exac_sev_hist::ExacerbationSeverityHistModule=ExacerbationSeverityHist(zeros(4), zeros(4)),
        total_hosp::Integer=0, family_hist::FamilyHistory=nothing, asthma_status::Bool=false,
        antibiotic_exposure::AntibioticExposureModule=nothing,
        province::String, month::Integer, SSP::String="SSP1_2.6"
    )

        census_division = CensusDivision(province, cal_year)
        pollution = Pollution(census_division.cduid, cal_year, month, SSP)

        if antibiotic_exposure !== nothing && family_hist !== nothing
            num_antibiotic_use = compute_num_antibiotic_use(
                antibiotic_exposure=antibiotic_exposure,
                sex=sex,
                birth_year=cal_year - age
            )
            has_family_hist = has_family_history_of_asthma(family_hist)
        end

        new(
            sex, age, cal_year, cal_year_index, alive, num_antibiotic_use, has_asthma,
            asthma_age, severity, control_levels, exac_hist, exac_sev_hist, total_hosp,
            has_family_hist, asthma_status, census_division, pollution
        )
    end
    function Agent(
        sex::Bool,
        age::Integer,
        cal_year::Integer,
        cal_year_index::Integer,
        alive::Bool,
        num_antibiotic_use::Integer,
        has_asthma::Bool,
        asthma_age::Union{Nothing, Integer},
        severity::Union{Nothing, Integer},
        control_levels::Union{Nothing, AbstractDict},
        exac_hist::Union{Nothing, ExacerbationHistModule},
        exac_sev_hist::Union{Nothing, ExacerbationSeverityHistModule},
        total_hosp::Integer,
        has_family_hist::Bool,
        asthma_status::Bool,
        census_division::Union{Nothing, CensusDivisionModule},
        pollution::Union{Nothing, PollutionModule}
    )
        new(
            sex, age, cal_year, cal_year_index, alive, num_antibiotic_use, has_asthma,
            asthma_age, severity, control_levels, exac_hist, exac_sev_hist, total_hosp,
            has_family_hist, asthma_status, census_division, pollution
        )
    end
end
