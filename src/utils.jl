function dict_initializer(parameter_names::Union{Nothing,Vector{Symbol}})
    isnothing(parameter_names) ? nothing : Dict(parameter_names .=> missing)
end

function vec_to_dict(v::AbstractArray, ll::AbstractVector)::AbstractDict
    d = Dict()
        for i in eachindex(ll)
            d[ll[i]] = v[i]
        end
    return d
end


function string_to_symbols_dict(dict::AbstractDict)::AbstractDict
    new_dict = Dict()
    for (key, value) in dict
        new_dict[Symbol(key)] = value
    end
    return new_dict
end


function set_up_birth(starting_year::Integer, population_growth_type::String, province::String)
    birth = Birth(nothing, nothing)
    @set! birth.estimate = filter(
        ([:year, :province, :projection_scenario] => (x, y, z) -> x >= starting_year
        && y == province && (z == population_growth_type || z == "past")),
        master_birth_estimate
    )
    relative(x) = x / birth.estimate.N[1]
    @set! birth.estimate = transform(birth.estimate, :N => relative)
    @set! birth.initial_population =  filter(
        ([:year, :province, :projection_scenario] => (x, y, z) -> x == starting_year
        && y == province && (z == population_growth_type || z == "past")),
        master_population_initial_distribution
    )
    return birth
end

function set_up_death(starting_year::Integer, province::String)
    death = Death(dict_initializer([:β0,:β1,:β2]), nothing)
    @set! death.parameters[:β0] = 0;
    @set! death.parameters[:β1] = 0;
    @set! death.parameters[:β2] = 0;
    @set! death.life_table = groupby(select(unstack(
        select(
            select(
                filter((
                    [:year, :province] => (x, y) -> x >= starting_year
                    && y == province
                    ),
                    master_life_table
                ),
            Not(:se)),
        Not(:province)),
    :sex,:prob_death),:F,:M,:year),:year)
    return death
end

function set_up_emigration(starting_year::Integer, population_growth_type::String,
    province::String)

    emigration = Emigration(nothing, nothing, nothing)
    @set! emigration.table = groupby(
        select(
            select(
                filter(
                    ([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year
                    && y == province
                    && z==population_growth_type),
                    master_emigration_table
                ),
                Not(:province)
            ),
            Not(:proj_scenario)
        ),
        :year
    )
    return emigration
end

function set_up_immigration(starting_year::Integer, population_growth_type::String,
    province::String)
    immigration = Immigration(nothing, nothing, nothing, nothing, nothing)
    @set! immigration.table = groupby(
        select(
            select(
                filter(
                    ([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year
                    && y == province
                    && z == population_growth_type),
                    master_immigration_table
                ),
                Not(:province)
            ),
            Not(:proj_scenario)
        ),
        :year
    )
    return immigration
end


function set_up_incidence(starting_year::Integer, province::String)::Incidence
    incidence = Incidence(
        dict_initializer([:β0_μ, :β0_σ]),
        dict_initializer([:β0]),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing
    )

    @set! incidence.hyperparameters[:β0_μ] = 0;
    @set! incidence.hyperparameters[:β0_σ] = 0.00000001;
    @set! incidence.incidence_table = groupby(
        filter([:year,:province] => (x,y) -> x >= min(
            starting_year, master_incidence_rate.year[nrow(master_incidence_rate)]
            ) && y==province, master_incidence_rate
            ),
        :year
    );
    @set! incidence.prevalence_table = groupby(
        filter([:year, :province] => (x,y) -> x >= min(
            starting_year-1, master_prevalence_rate.year[nrow(master_prevalence_rate)]
            ) && y==province, master_prevalence_rate),
        :year
    );
    @set! incidence.parameters[:β0] = 0;

    @set! incidence.calibration_table = groupby(
        select(
            filter([:province] => (x) -> x == province, M3_calibrated_asthma_prev_inc),
            Not([:province])
        ),
        [:year, :sex, :fam_history, :abx_exposure]
    );
    @set! incidence.min_year = collect(keys(incidence.calibration_table)[1])[1]+1
    @set! incidence.max_year = collect(
        keys(incidence.calibration_table)[length(incidence.calibration_table)]
        )[1]
    return incidence
end

function set_up_reassessment(starting_year::Integer, province::String)
    reassessment = Reassessment(nothing)
    @set! reassessment.table = groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_reassessment
        ),
        :year
    )
    return reassessment
end


function set_up_diagnosis(starting_year::Integer, province::String)
    diagnosis = Diagnosis(nothing, nothing)
    @set! diagnosis.true_positive_rates =  groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_dx
        ),
        :year
    )
    @set! diagnosis.false_negative_rates = groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_mis_dx
        ),
        :year
    )
    return diagnosis
end


function set_up_exacerbation_severity()
    exacerbation_severity = ExacerbationSeverity(
        dict_initializer([:α, :k]),
        dict_initializer([:p, :βprev_hosp_ped, :βprev_hosp_adult])
    )
    @set! exacerbation_severity.hyperparameters[:α] = [0.495, 0.195, 0.283, 0.026];
    @set! exacerbation_severity.hyperparameters[:k] = 100;
    @set! exacerbation_severity.parameters[:p] = ones(4) / 4;
    @set! exacerbation_severity.parameters[:βprev_hosp_ped] = 1.79
    @set! exacerbation_severity.parameters[:βprev_hosp_adult] = 2.88
    return exacerbation_severity
end


function set_up_utility()
    utility = Utility(dict_initializer([:eq5d, :control, :exac]))
    @set! utility.parameters[:eq5d] = eq5d
    # disutil
    @set! utility.parameters[:control] = [0.06, 0.09, 0.10]
    # disutil: duration 1 week for mild and two weeks for the rest
    @set! utility.parameters[:exac] = [0.32 * 1,  0.44 * 2, 0.50 * 2, 0.56 * 2 ] / 52
    return utility
end

function set_up_census_table()
    census_table = CensusTable(
        nothing, 2021
    )
    @set! census_table.data = groupby(master_census_data, :province)
    return census_table
end


function set_up(max_age=111, province="BC", starting_year=2000, time_horizon=19,
    num_births_initial=100, population_growth_type="LG")
    if province=="BC" || province=="CA"

        config = JSON.parsefile(CONFIG_PATH)

        agent = Agent(
            sex=false,
            age=0,
            cal_year=starting_year,
            cal_year_index=1,
            alive=true,
            num_antibiotic_use=0,
            has_asthma=false,
            asthma_age=0,
            severity=0,
            control_levels=nothing,
            exac_hist=ExacerbationHist(0, 0),
            exac_sev_hist=ExacerbationSeverityHist(zeros(4), zeros(4)),
            total_hosp=0,
            family_hist=false,
            asthma_status=false,
            census_division=nothing
        )

        birth = set_up_birth(starting_year, population_growth_type, province)
        death = set_up_death(starting_year, province)
        emigration = set_up_emigration(starting_year, population_growth_type, province)
        immigration = set_up_immigration(starting_year, population_growth_type, province)
        incidence = set_up_incidence(starting_year, province)
        reassessment = set_up_reassessment(starting_year, province)
        diagnosis = set_up_diagnosis(starting_year, province)
        exacerbation_severity = set_up_exacerbation_severity()
        utility = set_up_utility()
        census_table = set_up_census_table()

        simulation = Simulation(
            max_age=max_age,
            province=province,
            starting_calendar_year=starting_year,
            time_horizon=time_horizon,
            num_births_initial=num_births_initial,
            population_growth_type=population_growth_type,
            agent=agent,
            birth=birth,
            emigration=emigration,
            immigration=immigration,
            death=death,
            incidence=incidence,
            reassessment=reassessment,
            diagnosis=diagnosis,
            control=Control(config["control"]),
            exacerbation=Exacerbation(config["exacerbation"], province),
            exacerbation_severity=exacerbation_severity,
            antibiotic_exposure=AntibioticExposure(config["antibiotic_exposure"], abx_mid_trends, nothing),
            family_history=FamilyHistory(config["family_history"]),
            utility=utility,
            cost=AsthmaCost(config["cost"]),
            census_table=census_table,
            initial_distribution=nothing,
            outcome_matrix=(;)
        )

        return simulation
    else
        error("Province not supported")
    end
end
