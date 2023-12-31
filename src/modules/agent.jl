import Base.@kwdef


"""
    Agent

A person in the model.

# Fields
- `sex::Bool`: Sex of person, true = male, false = female.
- `age::Integer`: Age of person in years.
- `cal_year::Integer`: The calendar year, e.g. 2020.
- `cal_year_index::Integer`: TODO.
- `alive::Bool`: Whether the person is alive, true = alive.
- `num_antibiotic_use::Int`: TODO.
- `has_asthma::Bool`: Whether the person has astham, true = has asthma.
- `asthma_age::Union{Nothing, Int}`: Age at which the person was diagnosed with asthma.
- `severity::Union{Nothing, Int}`: Asthma severity level: 1 = mild, 2 = severe, 3 = very severe.
- `control::Union{Nothing,Vector{Float64}}`: Asthma control level: 1 = uncontrolled,
    2 = partially controlled, 3 = fully controlled.
- `exac_hist::Union{Nothing, ExacerbationHistModule}`: Total number of exacerbations.
- `exac_sev_hist::Union{Nothing, ExacerbationSeverityHistModule}`: Number of exacerbations by severity.
- `total_hosp::Integer`: Total number of very severe asthma exacerbations leading to hospitalization.
- `family_hist::Bool`: Is there a family history of asthma?
- `asthma_status::Bool`: TODO.
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
    control::Union{Nothing, Vector{Float64}}
    exac_hist::Union{Nothing, ExacerbationHistModule}
    exac_sev_hist::Union{Nothing, ExacerbationSeverityHistModule}
    total_hosp::Integer
    family_hist::Bool
    asthma_status::Bool
end

function set_agent!(ag, sex, age,cal_year,cal_year_index,alive, num_antibiotic_use,
    has_asthma,asthma_age,asthma_severity,asthma_control,asthma_exac_hist,asthma_exac_sev_history,total_hosp,fam_hist,asthma_status)
    ag.sex = sex
    ag.age= age
    ag.cal_year = cal_year
    ag.cal_year_index = cal_year_index
    ag.alive = alive
    ag.num_antibiotic_use=num_antibiotic_use
    ag.has_asthma = has_asthma
    ag.asthma_age = asthma_age
    ag.severity = asthma_severity
    ag.control = asthma_control
    ag.exac_hist = asthma_exac_hist
    ag.exac_sev_hist = asthma_exac_sev_history
    ag.total_hosp = total_hosp
    ag.family_hist = fam_hist
    ag.asthma_status = asthma_status
    nothing
end

function process_initial(ag::Agent,asthma_age_data)
    if ag.age==0
        return 0
    else
        return StatsBase.sample(Weights(asthma_age_data[1:ag.age+1,Int(ag.sex)+1]))-1
    end
end


"""
    create_agent(cal_year, cal_year_index, birth, age, antibiotic_exposure, family_hist, sex)

Creates a new agent (person).

# Arguments
- `cal_year::Integer`: the calendar year of the current iteration, e.g. 2027.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if
    the simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.
- `birth::Birth`: a Birth object, see  [`Birth`](@ref).
- `age::Integer`: the age of the person.
- `antibiotic_exposure::AntibioticExposure`: contains information about antibiotic exposure,
    see [`AntibioticExposure`](@ref).
- `family_hist::FamilyHistory`: contains information about family history of asthma,
    see [`FamilyHistory`](@ref).
- `sex::Bool`: sex of person, 1 = male, 0 = female, default = nothing.

# Returns
- `Agent`: a new agent.
"""
function create_agent(; cal_year::Integer, cal_year_index::Integer, sex::Bool, age::Integer,
    antibiotic_exposure::AntibioticExposureModule=nothing,
    family_hist::FamilyHistory=nothing)

    agent = Agent(
        sex=sex, age=age, cal_year=cal_year, cal_year_index=cal_year_index, alive=true,
        num_antibiotic_use=0, has_asthma=false, asthma_age=nothing, severity=nothing,
        control=nothing, exac_hist=ExacerbationHist(0, 0),
        exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)), total_hosp=0,
        family_hist=false, asthma_status=false
    )

    if antibiotic_exposure !== nothing && family_hist !== nothing
        if age == 0
            @set! agent.num_antibiotic_use = process_antibiotic_exposure(
                antibiotic_exposure, sex, cal_year
            )
            @set! agent.family_hist = process_family_history(
                family_hist
            )
        else
            @set! agent.num_antibiotic_use = process_antibiotic_exposure_initial(
                antibiotic_exposure, sex, cal_year - age
            )
            @set! agent.family_hist = process_family_history_initial(
                family_hist
            )
        end
    end
    return agent
end
