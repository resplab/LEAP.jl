"""
    Birth

A struct containing information about projected birth rates.

# Fields
- `estimate::Union{DataFrame, Nothing}`: A data frame giving the projected number of births in
    Canada with the following columns:
        `year`: integer year.
        `province`: A string indicating the province abbreviation, e.g. "BC".
        `N`: estimated number of births for that year.
        `prop_male`: proportion of births which are male, a number in [0, 1].
        `projection_scenario`: Population growth type, one of:
            ["past", "LG", "HG", "M1", "M2", "M3", "M4", "M5", "M6", FA", "SA"].
            See [Stats Canada](https://www150.statcan.gc.ca/n1/pub/91-520-x/91-520-x2022001-eng.htm).
        `N_relative`: number of births relative to the previous year.
    See `master_birth_estimate`.
- `initial_population::Union{DataFrame, Nothing}`: A data frame giving the population for the
    first year of the simulation:
        `year`: integer year.
        `age`: integer age.
        `province`: a string indicating the province abbreviation, e.g. "BC".
        `n`: estimated number of people in that age category in a given year.
        `n_birth`: the number of people born that year.
        `prop`: the ratio of that age group to the newborn age group (age = 0).
        `prop_male`: proportion of people in that age group who are male, a number in [0, 1].
        `projection_scenario`: Population growth type, one of:
            ["past", "LG", "HG", "M1", "M2", "M3", "M4", "M5", "M6", FA", "SA"].
            See [Stats Canada](https://www150.statcan.gc.ca/n1/pub/91-520-x/91-520-x2022001-eng.htm).
    See `master_population_initial_distribution`.
"""
struct Birth <: Birth_Module
    estimate::Union{DataFrame, Nothing}
    initial_population::Union{DataFrame, Nothing}
end

function create_agent_newborn(cal_year::Integer, cal_year_index::Integer, birth::Birth)
    agent = Agent(
        sex=rand(Bernoulli(birth.estimate.prop_male[cal_year_index])),
        age=0, cal_year=cal_year, cal_year_index=cal_year_index, alive=true, num_antibiotic_use=0,
        has_asthma=false, asthma_age=nothing, severity=nothing, control=nothing,
        exac_hist=ExacerbationHist(0, 0), exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)),
        total_hosp=0, family_hist=false, asthma_status=false
    )
    return agent
end


function create_agent_newborn(cal_year::Integer, cal_year_index::Integer, birth::Birth,
    antibiotic_exposure::AntibioticExposure, family_hist::FamilyHistory)
    agent = Agent(
        sex=rand(Bernoulli(birth.estimate.prop_male[cal_year_index])),
        age=0, cal_year=cal_year, cal_year_index=cal_year_index, alive=true, num_antibiotic_use=0,
        has_asthma=false, asthma_age=nothing, severity=nothing, control=nothing,
        exac_hist=ExacerbationHist(0, 0), exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)),
        total_hosp=0, family_hist=false, asthma_status=false
    )
    @set! agent.num_antibiotic_use = process(
        agent, antibiotic_exposure
    )
    @set! agent.family_hist = process(
        agent, family_hist
    )
    return agent
end


function create_agent_newborn(cal_year::Integer, cal_year_index::Integer, birth::Birth, sex::Bool,
    age::Integer, antibiotic_exposure::AntibioticExposure, family_hist::FamilyHistory)
    agent = Agent(
        sex=sex, age=age, cal_year=cal_year, cal_year_index=cal_year_index, alive=true,
        num_antibiotic_use=0, has_asthma=false, asthma_age=nothing, severity=nothing,
        control=nothing, exac_hist=ExacerbationHist(0, 0),
        exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)), total_hosp=0,
        family_hist=false, asthma_status=false
    )
    if age == 0
        @set! agent.num_antibiotic_use = process(
            agent, antibiotic_exposure
        )
        @set! agent.family_hist = process(
            agent, family_hist
        )
    else
        @set! agent.num_antibiotic_use = process_initial(
            agent,
            antibiotic_exposure,
            cal_year - agent.age
        )
        @set! agent.family_hist = process_initial(
            agent, family_hist
        )
    end
    return agent
end

function process_initial(b::Birth,n::Int)
    tmp_n = round.(Int,b.initial_population.prop*n)
    tmp_index = Vector{Int}[]
    for i in eachindex(tmp_n)
        tmp_index = vcat(tmp_index,fill(i,tmp_n[i]))
    end
    tmp_index
    # wsample(1:(nrow(b.initial_population)),b.initial_population.prop,n,replace=true)
end
