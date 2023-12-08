
mutable struct Simulation <: Simulation_Module
    max_age::Int
    province::Union{String,Char}
    starting_calendar_year::Int
    time_horizon::Union{Missing,Int,Vector{Int}}
    n::Union{Nothing,Missing,Real,String}
    population_growth_type::Union{Missing,String,Char}
    agent::Agent_Module
    birth::Birth_Module
    emigration::Emigration_Module
    immigration::Immigration_Module
    death::Death_Module
    incidence::Incidence_Module
    reassessment::Reassessment_Module
    diagnosis::Diagnosis_Module
    control::Control_Module
    exacerbation::Exacerbation_Module
    exacerbation_severity::Exacerbation_Severity_Module
    antibioticExposure::AntibioticExposure_Module
    familyHistory::FamilyHistory_Module
    util::Utility_Module
    cost::Cost_Module
    initial_distribution
    outcomeMatrix
end


function create_event_dict(until_all_die::Bool, cal_years::UnitRange{Int64}, min_cal_year::Integer,
    max_cal_year::Integer, max_age::Integer)
    event_list = [
        "antibiotic_exposure", "asthma_status", "asthma_incidence", "asthma_prevalence",
        "death","alive", "control", "exacerbation", "exacerbation_hospital",
        "exacerbation_by_severity", "emigration","immigration","family_history",
        "asthma_prevalence_family_history", "asthma_prevalence_antibiotic_exposure",
        "asthma_status_family_history", "asthma_status_antibiotic_exposure",
        "asthma_incidence_family_history", "asthma_incidence_antibiotic_exposure",
        "asthma_prevalence_contingency_table", "asthma_incidence_contingency_table",
        "family_history_prevalence", "util", "cost"
    ]

    event_dict = Dict()

    for event_name in event_list
        if event_name in ["asthma_prevalence_contingency_table", "asthma_incidence_contingency_table"]
            tmp_df = DataFrame(year=Int64[], sex=Int64[], age=Int64[], fam_history = Int64[],
                abx_exposure = Int64[], n_asthma= Int64[], n_no_asthma = Int64[])
            foreach(x -> push!(tmp_df, x), Iterators.product(
                min_cal_year:1:max_cal_year, 0:1:1, 0:1:max_age+1, 0:1:1, 0:1:3, 0, 0)
            )
            tmp_df = groupby(tmp_df, [:year, :sex, :fam_history, :abx_exposure])
            event_dict[event_name] = tmp_df
        else
            if event_name in ["control"]
                # year age sex 3-levels
                type = Real
                dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2, 3)
            elseif event_name in ["exacerbation_by_severity"]
                # year age sex 4-levels
                type = Real
                dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2, 4)
            elseif event_name in ["asthma_prevalence_family_history", "asthma_incidence_family_history",
                "asthma_status_family_history", "family_history_prevalence"]
                type = Int
                dimensions = (2, length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
            elseif event_name in ["asthma_prevalence_antibiotic_exposure",
                "asthma_incidence_antibiotic_exposure", "asthma_status_antibiotic_exposure"]
                type = Int
                dimensions = (4, length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
            elseif event_name in ["util","cost"]
                type = Real
                dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
            else
                # year age sex
                type = Int
                dimensions = (length(cal_years)+(until_all_die ? max_age : 0), max_age + 1, 2)
            end
            event_dict[event_name] = zeros(type, dimensions)
        end
    end
    return event_dict
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
    initial_pop_index = Int[]

    if cal_year == min_cal_year
        initial_pop_index = process_initial(simulation.birth, simulation.n)
        num_new_agents = length(initial_pop_index)
    end
    return num_new_agents
end


function generate_initial_asthma!(simulation::Simulation)
    @set! simulation.agent.has_asthma = agent_has_asthma(
        simulation.agent, simulation.incidence
    )
    if simulation.agent.has_asthma
        @set! simulation.agent.asthma_status = true
        @set! simulation.agent.asthma_age = compute_asthma_age(
            simulation.agent, simulation.incidence, simulation.agent.age
        )
        @set! simulation.agent.total_hosp = compute_hospitalization_prob(
            simulation.exacerbation_severity, simulation.agent.asthma_age, simulation
        )
        @set! simulation.agent.control = process_control(
            simulation.agent, simulation.control, true
        )
        @set! simulation.agent.exac_hist[1] = compute_num_exacerbations_initial(
            simulation.agent, simulation.exacerbation
        )
        # the number of exacerbation by severity
        @set! simulation.agent.exac_sev_hist[1] = process(
            simulation.exacerbation_severity, simulation.agent.exac_hist[1],
            (simulation.agent.total_hosp>0), simulation.agent.age
        )
        # update total hosp
        @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist[1][4]
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

    # if n is not provided, use a default value, n=100
    if ismissing(simulation.n)
        simulation.n = 100
    elseif simulation.n == "full"
        simulation.n = simulation.birth.initial_population.n_birth[1]
    elseif simulation.n < 1
        simulation.n = ceil(Int,simulation.n*simulation.birth.initial_population.n_birth[1])
    end

    max_time_horizon = (until_all_die ? typemax(Int) : simulation.time_horizon)

    # initiate some variables to collect results
    cal_years = min_cal_year:max_cal_year

    # store events
    n_list = zeros(Int,simulation.time_horizon,2)

    event_dict = create_event_dict(until_all_die, cal_years, min_cal_year, max_cal_year, max_age)

    # time the performance
    timer_output = TimerOutput()
    @timeit timer_output "sleep" sleep(0.02)

    # loop by year
    for cal_year in cal_years
        if verbose
            println(cal_year)
        end

        # time stamp
        @timeit timer_output "calendar year $cal_year" begin

        # index for cal_year
        tmp_cal_year_index = cal_year - min_cal_year + 1

        # num of newborns and immigrants in cal_year
        num_new_born = ceil(
            Int,
            simulation.n * simulation.birth.estimate.N_relative[tmp_cal_year_index]
        )
        num_immigrants = ceil(
            Int,
            num_new_born * sum(simulation.immigration.table[tmp_cal_year_index].n_prop_birth)
        )
        num_new_agents = get_num_new_agents(cal_year, min_cal_year, num_new_born, num_immigrants,
            simulation)

        # indicator for the new born;
        # otherwise immigrant
        new_born_indicator = vcat(falses(num_immigrants), trues(num_new_born))

        # weighted sampling of the immigrant profile
        if cal_year != min_cal_year
            immigrants_index = sample(1:nrow(simulation.immigration.table[tmp_cal_year_index]),
            Weights(simulation.immigration.table[tmp_cal_year_index].weights), num_immigrants)
        else
            initial_pop_index = process_initial(simulation.birth, simulation.n)
        end


        # for each agent i born/immigrated in cal_year
        for i in 1:num_new_agents

            # simulate an agent
            # random_parameter_initialization!(simulation.antibioticExposure)
            # random_parameter_initialization!(simulation.incidence)
            random_parameter_initialization!(simulation.control)
            random_parameter_initialization!(simulation.exacerbation)
            random_parameter_initialization!(simulation.exacerbation_severity)

            # if cal_year is the first year, generate the initial population
            if cal_year == min_cal_year
                tmp_index = initial_pop_index[i]
                simulation.agent = process_birth(cal_year, tmp_cal_year_index, simulation.birth,
                    rand(Bernoulli(simulation.birth.initial_population.prop_male[tmp_index])),
                    simulation.birth.initial_population.age[tmp_index],
                    simulation.antibioticExposure,
                    simulation.familyHistory
                )
            else # otherwise newborn or immigrant
                if new_born_indicator[i] # create a new-born
                    simulation.agent = process_birth(
                        cal_year, tmp_cal_year_index, simulation.birth,
                        simulation.antibioticExposure,
                        simulation.familyHistory
                    )
                else # immigrant
                    simulation.agent = process_immigration(
                        Bool(simulation.immigration.table[tmp_cal_year_index].sex[immigrants_index[i]]),
                        simulation.immigration.table[tmp_cal_year_index].age[immigrants_index[i]],
                        cal_year, tmp_cal_year_index, simulation.antibioticExposure,
                        simulation.familyHistory)
                    event_dict["immigration"][
                        simulation.agent.cal_year_index,
                        simulation.agent.age+1,
                        simulation.agent.sex+1] += 1
                end
            end

            event_dict["antibiotic_exposure"][simulation.agent.cal_year_index,
                simulation.agent.age + 1, simulation.agent.sex + 1] += simulation.agent.num_antibiotic_use
            event_dict["family_history"][simulation.agent.cal_year_index,
                simulation.agent.age + 1, simulation.agent.sex + 1] += simulation.agent.family_hist

            n_list[tmp_cal_year_index, simulation.agent.sex+1] +=1

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
                        simulation.agent, simulation.incidence, simulation.antibioticExposure.AbxOR
                    )
                    # crude incidence record
                    if simulation.agent.has_asthma
                        event_dict["asthma_incidence_contingency_table"][(
                            simulation.agent.cal_year, Int(simulation.agent.sex),
                            Int(simulation.agent.family_hist),
                            min(simulation.agent.num_antibiotic_use,3)
                            )][simulation.agent.age+1,"n_asthma"] += 1
                    else
                        event_dict["asthma_incidence_contingency_table"][(
                            simulation.agent.cal_year, Int(simulation.agent.sex),
                            Int(simulation.agent.family_hist),
                            min(simulation.agent.num_antibiotic_use,3)
                            )][simulation.agent.age+1,"n_no_asthma"] += 1
                    end

                    # keep track of patients who got asthma for the first time
                    if simulation.agent.has_asthma && !simulation.agent.asthma_status
                        @set! simulation.agent.asthma_status = true
                        @set! simulation.agent.asthma_age = simulation.agent.age
                        event_dict["asthma_status"][
                            simulation.agent.cal_year_index, simulation.agent.age+1,
                            simulation.agent.sex+1] += 1
                        # event_dict["asthma_status_family_history"][simulation.agent.family_hist+1,simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # event_dict["asthma_status_antibiotic_exposure"][min(simulation.agent.num_antibiotic_use+1,4),simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                    end

                    # asthma Dx
                    @set! simulation.agent.has_asthma = agent_has_asthma(
                        simulation.agent, simulation.diagnosis
                    )

                    # dx with asthma
                    if simulation.agent.has_asthma
                        # if they did not have asthma dx in the past, then record it
                        @set! simulation.agent.asthma_age = copy(simulation.agent.age)
                        event_dict["asthma_incidence"][
                            simulation.agent.cal_year_index, simulation.agent.age+1,
                            simulation.agent.sex+1] += 1
                        # event_dict["asthma_incidence_family_history"][simulation.agent.family_hist+1,simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # event_dict["asthma_incidence_antibiotic_exposure"][min(simulation.agent.num_antibiotic_use+1,4),simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1

                        # event_dict["asthma_prevalence"][simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # # event_dict["asthma_prevalence_family_history"][simulation.agent.family_hist+1,simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # # event_dict["asthma_prevalence_antibiotic_exposure"][min(simulation.agent.num_antibiotic_use+1,4),simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # event_dict["asthma_prevalence_contingency_table"][(simulation.agent.cal_year,Int(simulation.agent.sex),Int(simulation.agent.family_hist),min(simulation.agent.num_antibiotic_use,3))][simulation.agent.age+1,"n_asthma"] += 1

                        @set! simulation.agent.control = process(
                            simulation.agent, simulation.control
                        )
                        event_dict["control"][
                            simulation.agent.cal_year_index, simulation.agent.age+1,
                            simulation.agent.sex+1, :] += simulation.agent.control

                        @set! simulation.agent.exac_hist[1] = compute_num_exacerbations(
                            simulation.agent, simulation.exacerbation
                        )

                        if simulation.agent.exac_hist[1] != 0
                            @set! simulation.agent.exac_sev_hist[1] = process(
                                simulation.exacerbation_severity,
                                simulation.agent.exac_hist[1],
                                (simulation.agent.total_hosp>0),
                                simulation.agent.age
                            )
                            @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist[1][4]
                            event_dict["exacerbation"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1
                            ] += simulation.agent.exac_hist[1]
                            event_dict["exacerbation_hospital"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1
                            ] += simulation.agent.exac_sev_hist[1][4]
                            event_dict["exacerbation_by_severity"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1,:] .+= simulation.agent.exac_sev_hist[1]
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
                        # event_dict["asthma_prevalence"][simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # # event_dict["asthma_prevalence_family_history"][simulation.agent.family_hist+1,simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # # event_dict["asthma_prevalence_antibiotic_exposure"][min(simulation.agent.num_antibiotic_use+1,4),simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                        # event_dict["asthma_prevalence_contingency_table"][(simulation.agent.cal_year,Int(simulation.agent.sex),Int(simulation.agent.family_hist),min(simulation.agent.num_antibiotic_use,3))][simulation.agent.age+1,"n_asthma"] += 1

                        #  update control
                        @set! simulation.agent.control = process_control(
                            simulation.agent, simulation.control
                        )
                        event_dict["control"][
                            simulation.agent.cal_year_index,
                            simulation.agent.age+1,
                            simulation.agent.sex+1,:] += simulation.agent.control

                        # update exacerbation
                        @set! simulation.agent.exac_hist[2] = copy(simulation.agent.exac_hist[1])
                        @set! simulation.agent.exac_sev_hist[2] = copy(simulation.agent.exac_sev_hist[1])
                        @set! simulation.agent.exac_hist[1] = compute_num_exacerbations(
                            simulation.agent, simulation.exacerbation
                        )

                        if simulation.agent.exac_hist[1] != 0
                            @set! simulation.agent.exac_sev_hist[1] = process(
                                simulation.exacerbation_severity,
                                simulation.agent.exac_hist[1],
                                (simulation.agent.total_hosp>0),
                                simulation.agent.age
                            )
                            @set! simulation.agent.total_hosp += simulation.agent.exac_sev_hist[1][4]
                            event_dict["exacerbation"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1] += simulation.agent.exac_hist[1]
                            event_dict["exacerbation_hospital"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1] += simulation.agent.exac_sev_hist[1][4]
                            event_dict["exacerbation_by_severity"][
                                simulation.agent.cal_year_index,
                                simulation.agent.age+1,
                                simulation.agent.sex+1,:] .+= simulation.agent.exac_sev_hist[1]
                        end
                    end
                end

                # if no asthma, record it
                if simulation.agent.has_asthma
                    event_dict["asthma_prevalence"][simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += 1
                    event_dict["asthma_prevalence_contingency_table"][(simulation.agent.cal_year,Int(simulation.agent.sex),Int(simulation.agent.family_hist),min(simulation.agent.num_antibiotic_use,3))][simulation.agent.age+1,"n_asthma"] += 1
                else
                    event_dict["asthma_prevalence_contingency_table"][(simulation.agent.cal_year,Int(simulation.agent.sex),Int(simulation.agent.family_hist),min(simulation.agent.num_antibiotic_use,3))][simulation.agent.age+1,"n_no_asthma"] += 1
                end

                # util and cost
                event_dict["util"][simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += process(simulation.agent,simulation.util)
                event_dict["cost"][simulation.agent.cal_year_index,simulation.agent.age+1,simulation.agent.sex+1] += process(simulation.agent,simulation.cost)

                # death or emigration
                # assume death occurs first
                if process(simulation.agent,simulation.death)
                    @set! simulation.agent.alive = false
                    # everyone dies in the end... Inevitable
                    event_dict["death"][
                        simulation.agent.cal_year_index,
                        simulation.agent.age+1,
                        simulation.agent.sex+1
                    ] += 1
                # emigration
                elseif process(simulation.agent.cal_year_index,
                    simulation.agent.age,simulation.agent.sex,simulation.emigration)
                    @set! simulation.agent.alive = false
                    event_dict["emigration"][
                        simulation.agent.cal_year_index,
                        simulation.agent.age+1,
                        simulation.agent.sex+1] += 1
                else
                    # record alive
                    event_dict["alive"][
                        simulation.agent.cal_year_index,
                        simulation.agent.age+1,
                        simulation.agent.sex+1] += 1
                    # update the patient stats
                    @set! simulation.agent.age += 1
                    @set! simulation.agent.cal_year += 1
                    @set! simulation.agent.cal_year_index +=1
                end

            end # end while loop
        end # end for loop: agents
        # println(cal_year)
        end # end of begin timeit
    end # end for loop: cal year

    # return the output
    event_dict["asthma_prevalence_contingency_table"] = Matrix(combine(event_dict["asthma_prevalence_contingency_table"],[:year,:sex,:age,:fam_history,:abx_exposure,:n_asthma,:n_no_asthma]))
    event_dict["asthma_incidence_contingency_table"] = Matrix(combine(event_dict["asthma_incidence_contingency_table"],[:year,:sex,:age,:fam_history,:abx_exposure,:n_asthma,:n_no_asthma]))


    @set! simulation.outcomeMatrix = (; n = n_list, outcome_matrix = event_dict);

    if verbose
        print("\n Simulation finished. Check your simulation object for results.")
        print_timer(timer_output::TimerOutput)
    end

    return simulation.outcomeMatrix;
end
