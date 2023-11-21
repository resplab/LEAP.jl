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


"""
    set_agent!(agent, sex, age, cal_year, cal_year_index, alive, num_antibiotic_use,
        has_asthma, asthma_age, asthma_severity, asthma_control, asthma_exac_hist,
        asthma_exac_sev_history, total_hosp, fam_hist, asthma_status)

Given an agent object and parameters, initialize the values of the agent object.

# Arguments
- `agent::Agent`: An agent, see [`Agent`](@ref).
- `sex::Bool`: Sex of person, true = male, false = female.
- `age::Integer`: Age of person in years.
- `cal_year::Integer`: The calendar year, e.g. 2020.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.
- `alive::Bool`: Whether the person is alive, true = alive.
- `num_antibiotic_use::Integer`: TODO.
- `has_asthma::Bool`: Whether the person has astham, true = has asthma.
- `asthma_age::Union{Nothing, Integer}`: Age at which the person was diagnosed with asthma.
- `asthma_severity::Union{Nothing, Integer}`: Asthma severity level: 1 = mild, 2 = severe.
- `asthma_control::Union{Nothing, Vector{Integer}}`: Asthma control level: 1 = uncontrolled,
    2 = partially controlled, 3 = fully controlled.
- `asthma_exac_hist::Union{Nothing, Vector{Integer}}`: Total number of asthma exacerbations.
- `asthma_exac_sev_hist::Union{Nothing, Vector{Vector{Integer}}}`: Number of asthma exacerbations by
    severity.
- `total_hosp::Integer`: Total number of very severe asthma exacerbations leading to
    hospitalization.
- `fam_hist::Bool`: Is there a family history of asthma?
- `asthma_status::Bool`: TODO.

"""
function set_agent!(agent::Agent, sex::Bool, age::Integer, cal_year::Integer,
    cal_year_index::Integer, alive::Bool, num_antibiotic_use::Integer,
    has_asthma::Bool, asthma_age::Union{Nothing, Integer}, asthma_severity::Union{Nothing, Integer},
    asthma_control::Union{Nothing, Vector{Integer}},
    asthma_exac_hist::Union{Nothing, Vector{Integer}},
    asthma_exac_sev_history::Union{Nothing, Vector{Vector{Integer}}}, total_hosp::Integer,
    fam_hist::Bool, asthma_status::Bool)
    agent.sex = sex
    agent.age= age
    agent.cal_year = cal_year
    agent.cal_year_index = cal_year_index
    agent.alive = alive
    agent.num_antibiotic_use = num_antibiotic_use
    agent.has_asthma = has_asthma
    agent.asthma_age = asthma_age
    agent.severity = asthma_severity
    agent.control = asthma_control
    agent.exac_hist = asthma_exac_hist
    agent.exac_sev_hist = asthma_exac_sev_history
    agent.total_hosp = total_hosp
    agent.family_hist = fam_hist
    agent.asthma_status = asthma_status
    nothing
end


"""
    process_initial(agent, asthma_age_data)

TODO.

# Arguments
- `agent::Agent`: An agent, see [`Agent`](@ref).
- `asthma_age_data::Array{Float64, 2}`: An array of two vectors, with the first one corresponding
    to female and the second one to male. Each vector element represents the asthma incidence
    probability at the specified age.

"""
function process_initial(agent::Agent, asthma_age_data::Array{Float64, 2})
    if agent.age == 0
        return 0
    else
        return StatsBase.sample(Weights(asthma_age_data[1:agent.age + 1, Int(agent.sex) + 1])) - 1
    end
end
