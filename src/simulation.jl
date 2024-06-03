"""
    Simulation

TODO.

# Fields
- `max_age::Integer`: the maximum age to compute in the simulation.
- `province::Union{String, Char}`: a string indicating the province abbreviation, e.g. "BC".
- `min_year::Integer`: the calendar year to start the simulation at, e.g. 2000.
- `time_horizon::Union{Missing,Int,Vector{Int}}`: TODO.
- `num_births_initial::Union{Nothing,Missing,Real,String}`: the number of births for the initial
    year of the simulation.

"""
@kwdef mutable struct Simulation <: SimulationModule
    max_age::Integer
    province::Union{String,Char}
    min_year::Integer
    max_year::Integer
    time_horizon::Union{Missing,Int,Vector{Int}}
    num_births_initial::Union{Nothing,Missing,Real,String}
    population_growth_type::Union{Missing,String,Char}
    agent::Union{AgentModule, Nothing}
    birth::BirthModule
    emigration::EmigrationModule
    immigration::ImmigrationModule
    death::DeathModule
    incidence::IncidenceModule
    prevalence::PrevalenceModule
    reassessment::ReassessmentModule
    control::ControlModule
    exacerbation::ExacerbationModule
    exacerbation_severity::ExacerbationSeverityModule
    antibiotic_exposure::AntibioticExposureModule
    family_history::FamilyHistoryModule
    utility::UtilityModule
    cost::CostModule
    census_table::CensusTableModule
    pollution_table::PollutionTableModule
    SSP::String
    outcome_matrix
    function Simulation(config::AbstractDict)
        min_year = config["simulation"]["min_year"]
        max_age = config["simulation"]["max_age"]
        province = config["simulation"]["province"]
        population_growth_type = config["simulation"]["population_growth_type"]
        max_year = min_year + config["simulation"]["time_horizon"] - 1

        new(
            config["simulation"]["max_age"],
            province,
            min_year,
            max_year,
            config["simulation"]["time_horizon"],
            config["simulation"]["num_births_initial"],
            population_growth_type,
            nothing,
            Birth(min_year, province, population_growth_type, max_age),
            Emigration(min_year, province, population_growth_type),
            Immigration(min_year, province, population_growth_type, max_age),
            Death(config["death"], province, min_year),
            Incidence(config["incidence"]),
            Prevalence(config["prevalence"]),
            Reassessment(min_year, province),
            Control(config["control"]),
            Exacerbation(config["exacerbation"], province),
            ExacerbationSeverity(config["exacerbation_severity"]),
            AntibioticExposure(config["antibiotic_exposure"]),
            FamilyHistory(config["family_history"]),
            Utility(config["utility"]),
            AsthmaCost(config["cost"]),
            CensusTable(config["census_table"]),
            PollutionTable(),
            config["pollution"]["SSP"],
            (;)
        )
    end
    function Simulation(
        max_age::Integer,
        province::Union{String,Char},
        min_year::Integer,
        max_year::Integer,
        time_horizon::Union{Missing,Int,Vector{Int}},
        num_births_initial::Union{Nothing,Missing,Real,String},
        population_growth_type::Union{Missing,String,Char},
        agent::Union{AgentModule, Nothing},
        birth::BirthModule,
        emigration::EmigrationModule,
        immigration::ImmigrationModule,
        death::DeathModule,
        incidence::IncidenceModule,
        prevalence::PrevalenceModule,
        reassessment::ReassessmentModule,
        control::ControlModule,
        exacerbation::ExacerbationModule,
        exacerbation_severity::ExacerbationSeverityModule,
        antibiotic_exposure::AntibioticExposureModule,
        family_history::FamilyHistoryModule,
        utility::UtilityModule,
        cost::CostModule,
        census_table::CensusTableModule,
        pollution_table::PollutionTableModule,
        SSP::String,
        outcome_matrix
    )
        new(
            max_age,
            province,
            min_year,
            max_year,
            time_horizon,
            num_births_initial,
            population_growth_type,
            agent,
            birth,
            emigration,
            immigration,
            death,
            incidence,
            prevalence,
            reassessment,
            control,
            exacerbation,
            exacerbation_severity,
            antibiotic_exposure,
            family_history,
            utility,
            cost,
            census_table,
            pollution_table,
            SSP,
            outcome_matrix
        )
    end
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
    get_num_new_agents(year, min_year, num_new_born, num_immigrants, simulation)

TODO.

# Arguments
- `year::Integer`: the calendar year of the current iteration, e.g. 2027.
- `min_year::Integer`: the calendar year of the initial iteration, e.g. 2010.
- `num_new_born::Integer`: the number of babies born in the specified year `year`.
- `num_immigrants::Integer`: the number of immigrants who moved to Canada in the specified year
    `year`.
- `simulation::Simulation`:  Simulation module, see [`Simulation`](@ref).
"""
function get_num_new_agents(year::Integer, min_year::Integer, num_new_born::Integer,
    num_immigrants::Integer, simulation::Simulation)
    # for the first/initial year, we generate the initial population
    # otherwise we generate num_new_born + num_immigrants
    num_new_agents = (
        year==min_year ? ceil(
            Int, num_new_born / sum(filter(:age=> ==(0), simulation.birth.initial_population).prop)
        ) : num_new_born + num_immigrants
    )
    initial_pop_indices = Int[]

    if year == min_year
        initial_pop_indices = get_initial_population_indices(
            simulation.birth, simulation.num_births_initial
        )
        num_new_agents = length(initial_pop_indices)
    end
    return num_new_agents
end



function generate_initial_asthma!(simulation::Simulation)

    agent = deepcopy(simulation.agent)
    @set! agent.has_asthma = agent_has_asthma(agent, simulation.prevalence)

    if agent.has_asthma
        @set! agent.asthma_status = true
        @set! agent.asthma_age = compute_asthma_age(
            agent, simulation.incidence, simulation.prevalence, agent.age
        )
        @set! agent.total_hosp = compute_hospitalization_prob(
            agent, simulation.exacerbation_severity, simulation.control,
            simulation.exacerbation
        )
        @set! agent.control_levels = compute_control_levels(
            simulation.control, agent.sex, agent.age, true
        )
        @set! agent.exac_hist.num_current_year = compute_num_exacerbations_initial(
            agent, simulation.exacerbation
        )
        # the number of exacerbations by severity
        @set! agent.exac_sev_hist.current_year = compute_distribution_exac_severity(
            simulation.exacerbation_severity, agent.exac_hist.num_current_year,
            (agent.total_hosp>0), agent.age
        )
        # update total hospitalizations
        @set! agent.total_hosp += agent.exac_sev_hist.current_year[4]
    end
    setproperty!(simulation, Symbol("agent"), agent)
end


function get_new_agents(; simulation::SimulationModule, year::Integer, year_index::Integer)
    # num of newborns and immigrants in year
    num_new_born = get_num_newborn(
        simulation.birth, simulation.num_births_initial, year_index
    )
    num_immigrants = get_num_new_immigrants(
        simulation.immigration, num_new_born, year_index
    )
    num_new_agents = get_num_new_agents(year, simulation.min_year, num_new_born,
        num_immigrants, simulation
    )

    if year == simulation.min_year
        initial_pop_indices = get_initial_population_indices(
            simulation.birth, simulation.num_births_initial
        )
        sexes = [
            rand(
                Bernoulli(simulation.birth.initial_population.prop_male[index])
            )
            for index in initial_pop_indices
        ]
        ages = [
            simulation.birth.initial_population.age[index] for index in initial_pop_indices
        ]
        new_agents_df = DataFrame(age=ages, sex=sexes, immigrant=falses(num_new_agents))
    else
        immigrant_indices = sample(
            1:nrow(simulation.immigration.table[year_index]),
            Weights(simulation.immigration.table[year_index].prop_immigrants_year),
            num_immigrants
        )
        sexes_immigrant = [
            Bool(simulation.immigration.table[year_index].sex[index])
            for index in immigrant_indices
        ]
        sexes_birth = [
            rand(Bernoulli(simulation.birth.estimate.prop_male[year_index]))
            for index in (num_immigrants + 1):num_new_agents
        ]
        sexes = vcat(sexes_immigrant, sexes_birth)
        ages_immigrant = [
            simulation.immigration.table[year_index].age[index]
            for index in immigrant_indices
        ]
        ages_birth = [0 for index in (num_immigrants + 1):num_new_agents]
        ages = vcat(ages_immigrant, ages_birth)
        new_agents_df = DataFrame(
            age=ages, sex=sexes,
            immigrant=vcat(trues(num_immigrants), falses(num_new_born))
        )
    end
    return new_agents_df
end


function update_asthma_effects!(simulation::SimulationModule, outcome_matrix::OutcomeMatrixModule)

    agent = deepcopy(simulation.agent)
    @set! agent.control_levels = compute_control_levels(
        simulation.control, agent.sex, agent.age
    )
    add_control_to_outcome_matrix!(
        outcome_matrix,
        agent.age,
        agent.sex,
        agent.year_index,
        agent.control_levels
    )

    @set! agent.exac_hist.num_current_year = compute_num_exacerbations(
        agent, simulation.exacerbation
    )

    if agent.exac_hist.num_current_year != 0
        @set! agent.exac_sev_hist.current_year = (
        compute_distribution_exac_severity(
            simulation.exacerbation_severity,
            agent.exac_hist.num_current_year,
            (agent.total_hosp > 0),
            agent.age
        ))
        @set! agent.total_hosp += agent.exac_sev_hist.current_year[4]
        increment_field_in_outcome_matrix!(
            outcome_matrix,
            "exacerbation",
            agent.age,
            agent.sex,
            agent.year_index,
            agent.exac_hist.num_current_year
        )
        increment_field_in_outcome_matrix!(
            outcome_matrix,
            "exacerbation_hospital",
            agent.age, agent.sex, agent.year_index,
            agent.exac_sev_hist.current_year[4]
        )
        add_exacerbation_by_severity_to_outcome_matrix!(
            outcome_matrix,
            agent.age,
            agent.sex,
            agent.year_index,
            agent.exac_sev_hist.current_year
        )
    end
    setproperty!(simulation, Symbol("agent"), agent)
end


"""
    check_if_agent_gets_new_asthma_diagnosis!(simulation, outcome_matrix)

If the agent does not have asthma, check to see if they get a new diagnosis this year.
Mutates both the `simulation` and `outcome_matrix` arguments.

# Arguments
- `simulation::SimulationModule`: Simulation struct, see [`Simulation`](@ref).
- `outcome_matrix::OutcomeMatrixModule`: OutcomeMatrix struct, see [`OutcomeMatrix`](@ref).
"""
function check_if_agent_gets_new_asthma_diagnosis!(
    simulation::SimulationModule, outcome_matrix::OutcomeMatrixModule
)
    agent = deepcopy(simulation.agent)
    @set! agent.has_asthma = agent_has_asthma(
        agent, simulation.incidence, simulation.prevalence
    )
    # simulate and record asthma related events if they are labeled with asthma
    if agent.has_asthma
        # if they did not have asthma dx in the past, then record it
        @set! agent.asthma_age = copy(agent.age)
        increment_field_in_outcome_matrix!(outcome_matrix, "asthma_incidence",
            agent.age, agent.sex, agent.year_index
        )
        setproperty!(simulation, Symbol("agent"), agent)
        update_asthma_effects!(simulation, outcome_matrix)
        agent = deepcopy(simulation.agent)

        # keep track of patients who got asthma for the first time
        if !agent.asthma_status
            @set! agent.asthma_status = true
            increment_field_in_outcome_matrix!(outcome_matrix, "asthma_status",
                agent.age, agent.sex,
                agent.year_index
            )
        end
    end

    update_asthma_in_contingency_table!(outcome_matrix,
        agent.age, simulation.agent.sex,
        agent.year, agent.has_family_hist,
        agent.num_antibiotic_use,
        agent.has_asthma,
        "incidence"
    )
    setproperty!(simulation, Symbol("agent"), agent)
end


"""
    reassess_asthma_diagnosis!(simulation, outcome_matrix)

Reassess if the agent has asthma.

# Arguments
- `simulation::SimulationModule`: Simulation struct, see [`Simulation`](@ref).
- `outcome_matrix::OutcomeMatrixModule`: OutcomeMatrix struct, see [`OutcomeMatrix`](@ref).
"""
function reassess_asthma_diagnosis!(simulation::SimulationModule,
    outcome_matrix::OutcomeMatrixModule
)
    agent = deepcopy(simulation.agent)
    @set! agent.has_asthma = agent_has_asthma(
        agent, simulation.reassessment
    )
    # if still labeled with asthma
    if agent.has_asthma
        # update exacerbation
        @set! agent.exac_hist.num_prev_year = copy(
            agent.exac_hist.num_current_year
        )
        @set! agent.exac_sev_hist.prev_year = copy(
            agent.exac_sev_hist.current_year
        )
        setproperty!(simulation, Symbol("agent"), agent)
        update_asthma_effects!(simulation, outcome_matrix)
    else
        setproperty!(simulation, Symbol("agent"), agent)
    end
end



"""
    run_simulation(seed, until_all_die, verbose)

TODO.

# Arguments
- `seed:: Union{Missing, Float64}`: A seed to use for reproducibility.
- `until_all_die::Bool`: If true, run the simulation until all the agents have died. If false,
    use the `time_horizon` parameter and end the simulation after that number of years.
- `verbose::Bool`: If true, print out updates during simulation. Default true.
"""
function run_simulation(; seed=missing, until_all_die::Bool=false, verbose::Bool=false,
    config::Union{AbstractDict, Nothing}=nothing
)

    log_level = Logging.Warn
    if verbose
        log_level = Logging.Info
    end

    logger = ConsoleLogger(stdout, log_level)
    global_logger(logger)

    if isnothing(config)
        config = JSON.parsefile(CONFIG_PATH)
    end

    simulation = Simulation(config)

    if !ismissing(seed)
        Random.seed!(seed)
    end

    month = 1
    max_age = simulation.max_age
    min_year = simulation.min_year
    max_year = simulation.max_year

    max_time_horizon = (until_all_die ? typemax(Int) : simulation.time_horizon)
    years = min_year:max_year
    total_years = max_year - min_year + 1

    outcome_matrix = create_outcome_matrix(
        until_all_die=until_all_die,
        years=years,
        min_year=min_year,
        max_year=max_year,
        max_age=max_age
    )

    # time the performance
    timer_output = TimerOutput()
    @timeit timer_output "sleep" sleep(0.02)

    # loop by year
    for year in years
        @timeit timer_output "calendar year $year" begin

        year_index = year - min_year + 1

        @timeit timer_output "get_new_agents" begin
        new_agents_df = get_new_agents(
            simulation=simulation,
            year=year,
            year_index=year_index
        )
        end

        @info "Year $year, year $year_index of $total_years years, " *
        "$(size(new_agents_df)[1]) new agents born/immigrated."
        @info "$new_agents_df"
        # for each agent i born/immigrated in year
        for i in 1:size(new_agents_df)[1]
            assign_random_β0!(simulation.control)
            assign_random_β0!(simulation.exacerbation)
            assign_random_p!(simulation.exacerbation_severity)

            @timeit timer_output "create new Agent" begin
            simulation.agent = Agent(
                sex=new_agents_df.sex[i],
                age=new_agents_df.age[i],
                year=year,
                year_index=year_index,
                family_hist=simulation.family_history,
                antibiotic_exposure=simulation.antibiotic_exposure,
                province=simulation.province,
                month=month,
                SSP=simulation.SSP,
                census_table=simulation.census_table,
                pollution_table=simulation.pollution_table
            )
            end

            @info "[Agent #$(simulation.agent.uuid.short)]: $(i)th new agent of $year, " *
                  "age $(simulation.agent.age), sex $(Int(simulation.agent.sex)), " *
                  "immigrant: $(new_agents_df.immigrant[i]), " *
                  "newborn: $(!new_agents_df.immigrant[i])"

            @info "| -- Year: $(simulation.agent.year_index + min_year - 1), " *
            " age: $(simulation.agent.age)"

            if new_agents_df.immigrant[i]
                increment_field_in_outcome_matrix!(outcome_matrix, "immigration",
                    simulation.agent.age, simulation.agent.sex, simulation.agent.year_index
                )
            end

            increment_field_in_outcome_matrix!(outcome_matrix, "antibiotic_exposure",
                simulation.agent.age, simulation.agent.sex, simulation.agent.year_index,
                simulation.agent.num_antibiotic_use
            )

            increment_field_in_outcome_matrix!(outcome_matrix, "family_history",
                simulation.agent.age, simulation.agent.sex, simulation.agent.year_index,
                simulation.agent.has_family_hist
            )

            # if age >4, we need to generate the initial distribution of asthma related events
            if simulation.agent.age > 3
                @timeit timer_output "generate_initial_asthma!" begin
                generate_initial_asthma!(simulation)
                end
                @info "    | -- Agent age > 3, " *
                "agent has asthma (prevalence)? $(simulation.agent.has_asthma)"
            end

            # go through event processes for each agent
            while(simulation.agent.alive && simulation.agent.age <= max_age &&
                simulation.agent.year_index <= max_time_horizon)

                if !simulation.agent.has_asthma
                    @timeit timer_output "check_if_agent_gets_new_asthma_diagnosis!" begin
                    check_if_agent_gets_new_asthma_diagnosis!(simulation, outcome_matrix)
                    end
                    @info "    | -- Agent has asthma (incidence)? " *
                    "$(simulation.agent.has_asthma)"

                else
                    @timeit timer_output "reassess_asthma_diagnosis!" begin
                    reassess_asthma_diagnosis!(simulation, outcome_matrix)
                    end
                    @info "    | -- Agent was diagnosed with asthma, " *
                    "is this diagnosis correct? $(simulation.agent.has_asthma)"
                end

                # if no asthma, record it
                if simulation.agent.has_asthma
                    increment_field_in_outcome_matrix!(outcome_matrix, "asthma_prevalence",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.year_index
                    )
                end
                update_asthma_in_contingency_table!(outcome_matrix,
                    simulation.agent.age, simulation.agent.sex,
                    simulation.agent.year, simulation.agent.has_family_hist,
                    simulation.agent.num_antibiotic_use,
                    simulation.agent.has_asthma,
                    "prevalence"
                )

                @timeit timer_output "compute_utility" begin
                utility = compute_utility(simulation.agent, simulation.utility)
                @info "    | -- Utility of asthma: $utility"
                end

                increment_field_in_outcome_matrix!(
                    outcome_matrix, "utility", simulation.agent.age, simulation.agent.sex,
                    simulation.agent.year_index, utility
                )

                @timeit timer_output "compute_cost" begin
                cost = compute_cost(simulation.agent, simulation.cost)
                @info "    | -- Cost of asthma: $cost CAD"
                end

                increment_field_in_outcome_matrix!(
                    outcome_matrix, "cost", simulation.agent.age, simulation.agent.sex,
                    simulation.agent.year_index, cost
                )

                # death or emigration, assume death occurs first
                if agent_dies(simulation.agent, simulation.death)
                    @set! simulation.agent.alive = false
                    increment_field_in_outcome_matrix!(outcome_matrix, "death",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.year_index
                    )
                    @info "    | -- Agent has died at age $(simulation.agent.age)"
                # emigration
                elseif compute_prob_emigration(simulation.agent.year_index,
                    simulation.agent.age,simulation.agent.sex,simulation.emigration)
                    @set! simulation.agent.alive = false
                    increment_field_in_outcome_matrix!(outcome_matrix, "emigration",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.year_index
                    )
                    @info "    | -- Agent has emigrated at age $(simulation.agent.age)"
                else
                    # record alive
                    increment_field_in_outcome_matrix!(outcome_matrix, "alive",
                        simulation.agent.age, simulation.agent.sex, simulation.agent.year_index
                    )
                    # update the patient stats
                    @set! simulation.agent.age += 1
                    @set! simulation.agent.year += 1
                    @set! simulation.agent.year_index +=1

                    if (simulation.agent.age <= max_age &&
                        simulation.agent.year_index <= max_time_horizon
                    )
                        @info "| -- Year: $(simulation.agent.year_index + min_year - 1), " *
                        " age: $(simulation.agent.age)"
                    end
                end
            end
        end
        end
    end

    outcome_matrix.asthma_prevalence_contingency_table = Matrix(
            combine(outcome_matrix.asthma_prevalence_contingency_table,
            [:year, :sex, :age, :fam_history, :abx_exposure, :n_asthma, :n_no_asthma]
        )
    )
    outcome_matrix.asthma_incidence_contingency_table = Matrix(
        combine(
            outcome_matrix.asthma_incidence_contingency_table,
            [:year, :sex, :age, :fam_history, :abx_exposure, :n_asthma, :n_no_asthma]
        )
    )

    @set! simulation.outcome_matrix = outcome_matrix

    @info "\n Simulation finished. Check your simulation object for results."
    print_timer(timer_output::TimerOutput)

    return simulation.outcome_matrix
end
