function Dict_initializer(parameter_names::Union{Nothing,Vector{Symbol}})
    isnothing(parameter_names) ? nothing : Dict(parameter_names .=> missing)
end

function vec_to_dict(v::AbstractArray, ll::AbstractVector)::AbstractDict
    d = Dict()
        for i in eachindex(ll)
            d[ll[i]] = v[i]
        end
    return d
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
    death = Death(Dict_initializer([:β0,:β1,:β2]), nothing)
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
        Dict_initializer([:β0_μ, :β0_σ]),
        Dict_initializer([:β0]),
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


function set_up_diagnosis(starting_year::Integer, province::String)
    diagnosis = Diagnosis(nothing, nothing)
    @set! diagnosis.table =  groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_dx
        ),
        :year
    )
    @set! diagnosis.table_mis = groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_mis_dx
        ),
        :year
    )
    return diagnosis
end

function set_up_control()
    control = Control(
        Dict_initializer([:β0_μ, :β0_σ]),
        Dict_initializer([:β0, :βage, :βsex,:βsexage, :βsexage2, :βage2, :βDx2, :βDx3, :θ])
    )
    @set! control.hyperparameters[:β0_μ] = 0;
    @set! control.hyperparameters[:β0_σ] = 1.678728;
    @set! control.parameters[:βage] = 3.5430381;
    @set! control.parameters[:βage2] =-3.4980710;
    @set! control.parameters[:βsexage] = -0.8161495;
    @set! control.parameters[:βsexage2] = -1.1654264;
    @set! control.parameters[:βsex] = 0.2347807;
    @set! control.parameters[:θ] = [-0.3950; 2.754];
    return control
end


function set_up_exacerbation(province::String)
    exacerbation = Exacerbation(
        Dict_initializer([:β0_μ, :β0_σ]),
        Dict_initializer([:β0, :βage, :βsex, :βasthmaDx, :βprev_exac1, :βprev_exac2,
            :βcontrol_C, :βcontrol_PC, :βcontrol_UC, :calibration, :min_year]
        ),
        0
    )
    @set! exacerbation.initial_rate = 0.347;
    @set! exacerbation.hyperparameters[:β0_μ] = 0;
    @set! exacerbation.hyperparameters[:β0_σ] = 0.0000001;
    @set! exacerbation.parameters[:β0_calibration] = 0.0; # 0.056
    @set! exacerbation.parameters[:βage] = 0;
    @set! exacerbation.parameters[:βsex] = 0;
    @set! exacerbation.parameters[:βasthmaDx] = 0;
    @set! exacerbation.parameters[:βprev_exac1] = 0;
    @set! exacerbation.parameters[:βprev_exac2] = 0;
    @set! exacerbation.parameters[:βcontrol_C] =  log(0.1880058);
    @set! exacerbation.parameters[:βcontrol_PC] =  log(0.3760116);
    @set! exacerbation.parameters[:βcontrol_UC] =  log(0.5640174);
    @set! exacerbation.parameters[:βcontrol] =  0;
    @set! exacerbation.parameters[:calibration] = groupby(
        select(
            filter(
                [:province] => (x) -> x == province,
                exacerbation_calibration
            ),
            Not([:province])
        ),
        [:year,:sex]
    );
    @set! exacerbation.parameters[:min_year] = collect(
        keys(exacerbation.parameters[:calibration])[1]
    )[1]+1
    return exacerbation
end


function set_up_exacerbation_severity()
    exacerbation_severity = ExacerbationSeverity(
        Dict_initializer([:p0_μ, :p0_σ]),
        Dict_initializer([:p, :βprev_hosp_ped, :βprev_hosp_adult])
    )
    @set! exacerbation_severity.hyperparameters[:p0_μ] = [0.495, 0.195, 0.283, 0.026];
    @set! exacerbation_severity.hyperparameters[:p0_σ] = 100;
    @set! exacerbation_severity.parameters[:p] = ones(4) / 4;
    @set! exacerbation_severity.parameters[:βprev_hosp_ped] = 1.79
    @set! exacerbation_severity.parameters[:βprev_hosp_adult] = 2.88
    return exacerbation_severity
end


function set_up_antibiotic_exposure()
    antibiotic_exposure = AntibioticExposure(
        Dict_initializer([:β0_μ, :β0_σ]),
        Dict_initializer([:θ, :β0, :βage,:βsex, :βcal_year, :β2005, :β2005_cal_year,
            :fix2000, :βfloor, :midtrends]),
        nothing
    )
    @set! antibiotic_exposure.parameters[:θ] = 727.383;
    @set! antibiotic_exposure.parameters[:β0] = 110.000442;
    @set! antibiotic_exposure.parameters[:βage] = 0.0;
    @set! antibiotic_exposure.parameters[:βsex] = 0.249033;
    @set! antibiotic_exposure.parameters[:βcal_year] = -0.055100;
    @set! antibiotic_exposure.parameters[:β2005] = 55.033675;
    @set! antibiotic_exposure.parameters[:β2005_cal_year] = -0.027437;
    @set! antibiotic_exposure.parameters[:fixyear] = nothing;
    @set! antibiotic_exposure.parameters[:βfloor] = 50/1000;
    @set! antibiotic_exposure.parameters[:midtrends] = abx_mid_trends;
    return antibiotic_exposure
end


function set_up_family_history()
    family_history = FamilyHistory(nothing, Dict_initializer([:p]))
    @set! family_history.parameters[:p] = 0.2927242;
    return family_history
end


function set_up(max_age=111, province="BC", starting_year=2000, time_horizon=19, n=100,
    population_growth_type="LG")
    if province=="BC" || province=="CA"

        agent = Agent(
            false,
            0,
            starting_year,
            1,
            true,
            0,
            false,
            0,
            0,
            nothing,
            [0,0],
            [zeros(4), zeros(4)],
            0,
            false,
            false
        )

        birth = set_up_birth(starting_year, population_growth_type, province)
        death = set_up_death(starting_year, province)
        emigration = set_up_emigration(starting_year, population_growth_type, province)
        immigration = set_up_immigration(starting_year, population_growth_type, province)
        incidence = set_up_incidence(starting_year, province)

        reassessment = Reassessment(nothing)
        @set! reassessment.table =  groupby(filter([:year, :province] => (x, y) -> x >= starting_year && y == province, master_reassessment),:year)

        diagnosis = set_up_diagnosis(starting_year, province)
        control = set_up_control()
        exacerbation = set_up_exacerbation(province)
        exacerbation_severity = set_up_exacerbation_severity()
        antibiotic_exposure = set_up_antibiotic_exposure()
        family_history = set_up_family_history()

        util = Utility(Dict_initializer([:eq5d,:control,:exac]))
        @set! util.parameters[:eq5d] = eq5d
        # disutil
        @set! util.parameters[:control] = [0.06,0.09,0.10]
        # disutil: duration 1 week for mild and two weeks for the rest
        @set! util.parameters[:exac] = [0.32 * 1,  0.44 * 2 , 0.50 * 2 , 0.56 * 2 ] / 52

        cost = Cost(Dict_initializer([:control,:exac]))
        # 1.66 is the exchange rate btw 2018 USD and 2023 CAD Sept
        @set! cost.parameters[:control] = [2372, 2965, 3127]*1.66;
        @set! cost.parameters[:exac] = [130,594, 2425,9900]*1.66;

        simulation = Simulation(
            max_age,
            province,
            starting_year,
            time_horizon,
            n,
            population_growth_type,
            agent,
            birth,
            emigration,
            immigration,
            death,
            incidence,
            reassessment,
            diagnosis,
            control,
            exacerbation,
            exacerbation_severity,
            antibiotic_exposure,
            family_history,
            util,
            cost,
            nothing,
            (;)
        )

        return simulation
    else
        error("Province not supported")
    end
end
