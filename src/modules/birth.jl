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


"""
    get_initial_population_indices(birth, num_births)

Get the indices for the agents from the initial population table, weighted by age.

# Examples
For example, if the number of births is 2, and we have the following initial population table:

    age | prop | ...
    ----------------
    0     1.0    ...
    1     2.0    ...
    2     0.5    ...

then we will return the following:

    [1, 1, 2, 2, 2, 2, 3]

# Arguments
- `birth::Birth`: a Birth object, see  [`Birth`](@ref).
- `num_births::Integer`: number of births.

# Returns
- `Vector{Integer}`: the indices for the initial population table.
"""
function get_initial_population_indices(birth::Birth, num_births::Integer)
    num_agents_per_age_group = round.(Int, birth.initial_population.prop*num_births)
    intial_population_indices = Vector{Integer}[]
    for age_index in eachindex(num_agents_per_age_group)
        initial_population_indices = vcat(
            initial_population_indices,
            fill(age_index, num_agents_per_age_group[age_index])
        )
    end
    return intial_population_indices
end


"""
    get_num_newborn(birth, num_births_initial, cal_year_index)

Get the number of births in a given year.

# Arguments
- `birth::Birth`: a Birth object, see  [`Birth`](@ref).
- `num_births_initial::Integer`: number of births in the initial year of the simulation.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.

# Returns
- `Integer`: the number of births for the given year.
"""
function get_num_newborn(birth::Birth, num_births_initial::Integer, cal_year_index::Integer)
    num_new_born = ceil(
        Int,
        num_births_initial * birth.estimate.N_relative[cal_year_index]
    )
    return num_new_born
end
