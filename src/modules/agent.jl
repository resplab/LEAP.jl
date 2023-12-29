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
- `exac_hist::Union{Nothing, ExacerbationHist_Module}`: Total number of exacerbations.
- `exac_sev_hist::Union{Nothing, ExacerbationSeverityHist_Module}`: Number of exacerbations by severity.
- `total_hosp::Integer`: Total number of very severe asthma exacerbations leading to hospitalization.
- `family_hist::Bool`: Is there a family history of asthma?
- `asthma_status::Bool`: TODO.
- `census_division::CensusDivisionModule`: the Canadian census division where the agent resides.
"""
@kwdef struct Agent  <: Agent_Module
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
    exac_hist::Union{Nothing, ExacerbationHist_Module}
    exac_sev_hist::Union{Nothing, ExacerbationSeverityHist_Module}
    total_hosp::Integer
    family_hist::Bool
    asthma_status::Bool
    census_division::CensusDivisionModule
end

function set_agent!(agent, sex, age,cal_year,cal_year_index,alive, num_antibiotic_use,
    has_asthma,asthma_age,asthma_severity,asthma_control,asthma_exac_hist,asthma_exac_sev_history,
    total_hosp,fam_hist, asthma_status, census_division)
    agent.sex = sex
    agent.age= age
    agent.cal_year = cal_year
    agent.cal_year_index = cal_year_index
    agent.alive = alive
    agent.num_antibiotic_use=num_antibiotic_use
    agent.has_asthma = has_asthma
    agent.asthma_age = asthma_age
    agent.severity = asthma_severity
    agent.control = asthma_control
    agent.exac_hist = asthma_exac_hist
    agent.exac_sev_hist = asthma_exac_sev_history
    agent.total_hosp = total_hosp
    agent.family_hist = fam_hist
    agent.asthma_status = asthma_status
    agent.census_division = census_division
end

function process_initial(agent::Agent,asthma_age_data)
    if agent.age==0
        return 0
    else
        return StatsBase.sample(Weights(asthma_age_data[1:agent.age+1,Int(agent.sex)+1]))-1
    end
end
