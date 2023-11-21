"""
    Agent

A person in the model.

# Fields
- `sex::Bool`: Sex of person, true = male, false = female.
- `age::Integer`: Age of person in years.
- `cal_year::Integer`: The calendar year, e.g. 2020.
- `cal_year_index::Integer`: TODO.
- `alive::Bool`: Whether the person is alive, true = alive.
- `num_antibiotic_use::Integer`: TODO.
- `has_asthma::Bool`: Whether the person has astham, true = has asthma.
- `asthma_age::Union{Nothing, Integer}`: Age at which the person was diagnosed with asthma.
- `severity::Union{Nothing, Integer}`: Asthma severity level: 1 = mild, 2 = severe, 3 = very severe.
- `control::Union{Nothing,Vector{Float64}}`: Asthma control level: 1 = uncontrolled,
    2 = partially controlled, 3 = fully controlled.
- `exac_hist::Union{Nothing, Vector{Integer}}`: Total number of asthma exacerbations.
- `exac_sev_hist::Union{Nothing,Vector{Vector{Integer}}}`: Number of asthma exacerbations by severity.
- `total_hosp::Integer`: Total number of very severe asthma exacerbations leading to hospitalization.
- `family_hist::Bool`: Is there a family history of asthma?
- `asthma_status::Bool`: TODO.
"""
struct Agent  <: Agent_Module
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
    exac_hist::Union{Nothing, Vector{Integer}}
    exac_sev_hist::Union{Nothing, Vector{Vector{Integer}}}
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
