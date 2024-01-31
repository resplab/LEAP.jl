import Base.@kwdef

"""
    Simulation

TODO.

# Fields
- `max_age::Integer`: the maximum age to compute in the simulation.
- `province::Union{String, Char}`: a string indicating the province abbreviation, e.g. "BC".
- `starting_calendar_year::Integer`: the calendar year to start the simulation at, e.g. 2000.
- `time_horizon::Union{Missing,Int,Vector{Int}}`: TODO.
- `num_births_initial::Union{Nothing,Missing,Real,String}`: the number of births for the initial
    year of the simulation.

"""
@kwdef mutable struct Simulation <: SimulationModule
    max_age::Integer
    province::Union{String,Char}
    starting_calendar_year::Integer
    time_horizon::Union{Missing,Int,Vector{Int}}
    num_births_initial::Union{Nothing,Missing,Real,String}
    population_growth_type::Union{Missing,String,Char}
    agent::AgentModule
    birth::BirthModule
    emigration::EmigrationModule
    immigration::ImmigrationModule
    death::DeathModule
    incidence::IncidenceModule
    reassessment::ReassessmentModule
    diagnosis::DiagnosisModule
    control::ControlModule
    exacerbation::ExacerbationModule
    exacerbation_severity::ExacerbationSeverityModule
    antibioticExposure::AntibioticExposureModule
    familyHistory::FamilyHistoryModule
    util::UtilityModule
    cost::CostModule
    census_table::CensusTableModule
    initial_distribution
    outcome_matrix
end

function set_num_births_initial!(simulation::SimulationModule, num_births_initial::Integer=100)
    if num_births_initial == "full"
        simulation.num_births_initial = simulation.birth.initial_population.n_birth[1]
    elseif num_births_initial < 1
        simulation.num_births_initial = ceil(
            Int,
            num_births_initial*simulation.birth.initial_population.n_birth[1]
        )
    end
end


"""
    get_num_new_agents(cal_year, min_cal_year, num_new_born, num_immigrants, simulation)

TODO.

# Arguments
- `cal_year::Integer`: the calendar year of the current iteration, e.g. 2027.
- `min_cal_year::Integer`: the calendar year of the initial iteration, e.g. 2010.
- `num_new_born::Integer`: the number of babies born in the specified year `cal_year`.
- `num_immigrants::Integer`: the number of immigrants who moved to Canada in the specified year
    `cal_year`.
- `simulation::Simulation`:  Simulation module, see [`Simulation`](@ref).
"""
function get_num_new_agents(cal_year::Integer, min_cal_year::Integer, num_new_born::Integer,
    num_immigrants::Integer, simulation::Simulation)
    # for the first/initial year, we generate the initial population
    # otherwise we generate num_new_born + num_immigrants
    num_new_agents = (
        cal_year==min_cal_year ? ceil(
            Int, num_new_born / sum(filter(:age=> ==(0), simulation.birth.initial_population).prop)
        ) : num_new_born + num_immigrants
    )
    initial_pop_indices = Int[]

    if cal_year == min_cal_year
        initial_pop_indices = get_initial_population_indices(
            simulation.birth, simulation.num_births_initial
        )
        num_new_agents = length(initial_pop_indices)
    end
    return num_new_agents
end



function generate_initial_asthma!(simulation::Simulation)
    @set! simulation.agent.has_asthma = agent_has_asthma(
        simulation.agent, simulation.incidence, "prevalence"
    )
    if simulation.agent.has_asthma
        @set! simulation.agent.asthma_status = true
        @set! simulation.agent.asthma_age = compute_asthma_age(
            simulation.agent, simulation.incidence, simulation.agent.age
        )
        @set! simulation.agent.total_hosp = compute_hospitalization_prob(
            simulation.exacerbation_severity, simulation.agent.asthma_age, simulation
        )
        @set! simulation.agent.control = compute_control_levels_prob(
            simulation.control, simulation.agent.sex, simulation.agent.age, true
        )
        @set! simulation.agent.exac_hist.num_current_year = compute_num_exacerbations_initial(
            simulation.agent, simulation.exacerbation
        )
        # the number of exacerbation by severity
        @set! simulation.agent.exac_sev_hist.current_year = compute_distribution_exac_severity(
            simulation.exacerbation_severity, simulation.agent.exac_hist.num_current_year,
            (simulation.agent.total_hosp>0), simulation.agent.age
        )
        # update total hosp
        @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist.current_year[4]
    end
end


"""
    process(simulation, seed, until_all_die, verbose)

TODO.

# Arguments
- `simulation::Simulation`:  Simulation module, see [`Simulation`](@ref).
- `seed:: Union{Missing, Float64}`: TODO.
- `until_all_die::Bool`: TODO.
- `verbose::Bool`: If true, print out updates during simulation. Default true.
"""
function process(simulation::Simulation, seed=missing, until_all_die::Bool=false, verbose::Bool=true)
    # reproducibility
    if !ismissing(seed)
        Random.seed!(seed)
    end

    max_age = simulation.max_age
    min_cal_year = simulation.starting_calendar_year
    max_cal_year = min_cal_year + simulation.time_horizon - 1

    max_time_horizon = (until_all_die ? typemax(Int) : simulation.time_horizon)
    cal_years = min_cal_year:max_cal_year

    # store events
    n_list = zeros(Int, simulation.time_horizon, 2)

    outcome_matrix = create_outcome_matrix(until_all_die, cal_years, min_cal_year, max_cal_year, max_age)

    # time the performance
    timer_output = TimerOutput()
    @timeit timer_output "sleep" sleep(0.02)

    # loop by year
    for cal_year in cal_years
        # time stamp
        @timeit timer_output "calendar year $cal_year" begin

        if verbose
            println(cal_year)
        end

        # index for cal_year
        tmp_cal_year_index = cal_year - min_cal_year + 1

        # num of newborns and immigrants in cal_year
        num_new_born = get_num_newborn(
            simulation.birth, simulation.num_births_initial, tmp_cal_year_index
        )
        num_immigrants = get_num_new_immigrants(
            simulation.immigration, num_new_born, tmp_cal_year_index
        )
        num_new_agents = get_num_new_agents(cal_year, min_cal_year, num_new_born, num_immigrants,
            simulation)

        # indicator for the new born; otherwise immigrant
        new_born_indicator = vcat(falses(num_immigrants), trues(num_new_born))

        # weighted sampling of the immigrant profile
        if cal_year != min_cal_year
            immigrant_indices = sample(
                1:nrow(simulation.immigration.table[tmp_cal_year_index]),
                Weights(simulation.immigration.table[tmp_cal_year_index].weights),
                num_immigrants
            )
        else
            initial_pop_indices = get_initial_population_indices(
                simulation.birth, simulation.num_births_initial
            )
        end


        # for each agent i born/immigrated in cal_year
        for i in 1:num_new_agents

            # simulate an agent
            random_parameter_initialization!(simulation.control)
            random_parameter_initialization!(simulation.exacerbation)
            random_parameter_initialization!(simulation.exacerbation_severity)

            if cal_year == min_cal_year
                sex = rand(
                    Bernoulli(simulation.birth.initial_population.prop_male[initial_pop_indices[i]])
                )
                age = simulation.birth.initial_population.age[initial_pop_indices[i]]
                simulation.agent = create_agent(
                    cal_year=cal_year,
                    cal_year_index=tmp_cal_year_index,
                    sex=sex,
                    age=age,
                    province=simulation.province,
                    antibiotic_exposure=simulation.antibioticExposure,
                    family_hist=simulation.familyHistory,
                    census_table=simulation.census_table
                )
            elseif new_born_indicator[i]
                simulation.agent = create_agent(
                    cal_year=cal_year,
                    cal_year_index=tmp_cal_year_index,
                    sex=rand(Bernoulli(simulation.birth.estimate.prop_male[tmp_cal_year_index])),
                    age=0,
                    province=simulation.province,
                    antibiotic_exposure=simulation.antibioticExposure,
                    family_hist=simulation.familyHistory,
                    census_table=simulation.census_table
                )
            else
                simulation.agent = create_agent(
                    cal_year=cal_year,
                    cal_year_index=tmp_cal_year_index,
                    sex=Bool(simulation.immigration.table[tmp_cal_year_index].sex[immigrant_indices[i]]),
                    age=simulation.immigration.table[tmp_cal_year_index].age[immigrant_indices[i]],
                    province=simulation.province,
                    antibiotic_exposure=simulation.antibioticExposure,
                    family_hist=simulation.familyHistory,
                    census_table=simulation.census_table
                )
                increment_field_in_outcome_matrix!(outcome_matrix, "immigration",
                    simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                )
            end

            increment_field_in_outcome_matrix!(outcome_matrix, "antibiotic_exposure",
                simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index,
                simulation.agent.family_hist
            )

            increment_field_in_outcome_matrix!(outcome_matrix, "family_history",
                simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index,
                simulation.agent.num_antibiotic_use
            )

            n_list[tmp_cal_year_index, simulation.agent.sex+1] += 1

            # if age >4, we need to generate the initial distribution of asthma related events
            if simulation.agent.age > 3
                generate_initial_asthma!(simulation)
            end

            # go through event processes for each agent
            while(simulation.agent.alive && simulation.agent.age <= max_age &&
                simulation.agent.cal_year_index <= max_time_horizon)
                # no asthma
                if !simulation.agent.has_asthma
                    # asthma inc
                    @set! simulation.agent.has_asthma = agent_has_asthma(
                        simulation.agent, simulation.incidence, "incidence"
                    )
                    # crude incidence record
                    if simulation.agent.has_asthma
                        # keep track of patients who got asthma for the first time
                        if !simulation.agent.asthma_status
                            @set! simulation.agent.asthma_status = true
                            @set! simulation.agent.asthma_age = simulation.agent.age
                            increment_field_in_outcome_matrix!(outcome_matrix, "asthma_status",
                                simulation.agent.age, simulation.agent.sex,
                                simulation.agent.cal_year_index
                            )
                        end
                    end

                    update_asthma_in_contingency_table!(outcome_matrix,
                        simulation.agent.age, simulation.agent.sex,
                        simulation.agent.cal_year, simulation.agent.family_hist,
                        simulation.agent.num_antibiotic_use,
                        simulation.agent.has_asthma,
                        "incidence"
                    )

                    # asthma Dx
                    @set! simulation.agent.has_asthma = agent_has_asthma(
                        simulation.agent, simulation.diagnosis
                    )

                    # dx with asthma
                    if simulation.agent.has_asthma
                        # if they did not have asthma dx in the past, then record it
                        @set! simulation.agent.asthma_age = copy(simulation.agent.age)
                        increment_field_in_outcome_matrix!(outcome_matrix, "asthma_incidence",
                            simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                        )

                        @set! simulation.agent.control = compute_control_levels_prob(
                            simulation.control, simulation.agent.sex, simulation.agent.age
                        )
                        add_control_to_outcome_matrix!(outcome_matrix, simulation.agent.age,
                            simulation.agent.sex, simulation.agent.cal_year_index,
                            simulation.agent.control
                        )

                        @set! simulation.agent.exac_hist.num_current_year = compute_num_exacerbations(
                            simulation.agent, simulation.exacerbation
                        )

                        if simulation.agent.exac_hist.num_current_year != 0
                            @set! simulation.agent.exac_sev_hist.current_year = (
                            compute_distribution_exac_severity(
                                simulation.exacerbation_severity,
                                simulation.agent.exac_hist.num_current_year,
                                (simulation.agent.total_hosp>0),
                                simulation.agent.age
                            ))
                            @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist.current_year[4]
                            increment_field_in_outcome_matrix!(
                                outcome_matrix,
                                "exacerbation",
                                simulation.agent.age,
                                simulation.agent.sex,
                                simulation.agent.cal_year_index
                            )
                            increment_field_in_outcome_matrix!(outcome_matrix, "exacerbation_hospital",
                                simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index,
                                simulation.agent.exac_sev_hist.current_year[4]
                            )
                            add_exacerbation_by_severity_to_outcome_matrix!(outcome_matrix,
                                simulation.agent.age,
                                simulation.agent.sex, simulation.agent.cal_year_index,
                                simulation.agent.exac_sev_hist.current_year
                            )
                        end
                    end
                # has asthma
                else
                    # reassessment
                    @set! simulation.agent.has_asthma = agent_has_asthma(
                        simulation.agent, simulation.reassessment
                    )
                    # if still dxed with asthma
                    if simulation.agent.has_asthma
                        #  update control
                        @set! simulation.agent.control = compute_control_levels_prob(
                            simulation.control, simulation.agent.sex, simulation.agent.age
                        )
                        add_control_to_outcome_matrix!(outcome_matrix, simulation.agent.age,
                            simulation.agent.sex, simulation.agent.cal_year_index,
                            simulation.agent.control
                        )

                        # update exacerbation
                        @set! simulation.agent.exac_hist.num_prev_year = copy(
                            simulation.agent.exac_hist.num_current_year
                        )
                        @set! simulation.agent.exac_sev_hist.prev_year = copy(
                            simulation.agent.exac_sev_hist.current_year
                        )
                        @set! simulation.agent.exac_hist.num_current_year = compute_num_exacerbations(
                            simulation.agent, simulation.exacerbation
                        )

                        if simulation.agent.exac_hist.num_current_year != 0
                            @set! simulation.agent.exac_sev_hist.current_year = (
                            compute_distribution_exac_severity(
                                simulation.exacerbation_severity,
                                simulation.agent.exac_hist.num_current_year,
                                (simulation.agent.total_hosp>0),
                                simulation.agent.age
                            ))
                            @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist.current_year[4]
                            increment_field_in_outcome_matrix!(
                                outcome_matrix,
                                "exacerbation",
                                simulation.agent.age,
                                simulation.agent.sex,
                                simulation.agent.cal_year_index
                            )
                            increment_field_in_outcome_matrix!(
                                outcome_matrix,
                                "exacerbation_hospital",
                                simulation.agent.age,
                                simulation.agent.sex,
                                simulation.agent.cal_year_index,
                                simulation.agent.exac_sev_hist.current_year[4]
                            )
                            add_exacerbation_by_severity_to_outcome_matrix!(outcome_matrix,
                                simulation.agent.age,
                                simulation.agent.sex, simulation.agent.cal_year_index,
                                simulation.agent.exac_sev_hist.current_year
                            )
                        end
                    end
                end

                # if no asthma, record it
                if simulation.agent.has_asthma
                    increment_field_in_outcome_matrix!(outcome_matrix, "asthma_prevalence",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                    )
                end
                update_asthma_in_contingency_table!(outcome_matrix,
                    simulation.agent.age, simulation.agent.sex,
                    simulation.agent.cal_year, simulation.agent.family_hist,
                    simulation.agent.num_antibiotic_use,
                    simulation.agent.has_asthma,
                    "prevalence"
                )

                # util and cost
                increment_field_in_outcome_matrix!(outcome_matrix, "util",
                    simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index,
                    compute_utility(simulation.agent, simulation.util)
                )
                increment_field_in_outcome_matrix!(outcome_matrix, "cost",
                    simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index,
                    compute_cost(simulation.agent, simulation.cost)
                )

                # death or emigration, assume death occurs first
                if compute_prob_death(simulation.agent, simulation.death)
                    @set! simulation.agent.alive = false
                    increment_field_in_outcome_matrix!(outcome_matrix, "death",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                    )
                # emigration
                elseif compute_prob_emigration(simulation.agent.cal_year_index,
                    simulation.agent.age,simulation.agent.sex,simulation.emigration)
                    @set! simulation.agent.alive = false
                    increment_field_in_outcome_matrix!(outcome_matrix, "emigration",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                    )
                else
                    # record alive
                    increment_field_in_outcome_matrix!(outcome_matrix, "alive",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.cal_year_index
                    )
                    # update the patient stats
                    @set! simulation.agent.age += 1
                    @set! simulation.agent.cal_year += 1
                    @set! simulation.agent.cal_year_index +=1
                end

            end
        end
        end
    end

    outcome_matrix.asthma_prevalence_contingency_table = Matrix(
            combine(outcome_matrix.asthma_prevalence_contingency_table,
            [:year,:sex,:age,:fam_history,:abx_exposure,:n_asthma,:n_no_asthma]
        )
    )
    outcome_matrix.asthma_incidence_contingency_table = Matrix(
        combine(
            outcome_matrix.asthma_incidence_contingency_table,
            [:year,:sex,:age,:fam_history,:abx_exposure,:n_asthma,:n_no_asthma]
        )
    )


    @set! simulation.outcome_matrix = (; n = n_list, outcome_matrix = outcome_matrix);

    if verbose
        print("\n Simulation finished. Check your simulation object for results.")
        print_timer(timer_output::TimerOutput)
    end

    return simulation.outcome_matrix;
end
