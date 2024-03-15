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

        reassessment = set_up_reassessment(starting_year, province)
        diagnosis = set_up_diagnosis(starting_year, province)

        simulation = Simulation(
            max_age=max_age,
            province=province,
            starting_calendar_year=starting_year,
            time_horizon=time_horizon,
            num_births_initial=num_births_initial,
            population_growth_type=population_growth_type,
            agent=agent,
            birth=Birth(starting_year, province, population_growth_type),
            emigration=Emigration(starting_year, province, population_growth_type),
            immigration=Immigration(starting_year, province, population_growth_type),
            death=Death(config["death"], province, starting_year),
            incidence=Incidence(config["incidence"], starting_year, province),
            reassessment=reassessment,
            diagnosis=diagnosis,
            control=Control(config["control"]),
            exacerbation=Exacerbation(config["exacerbation"], province),
            exacerbation_severity=ExacerbationSeverity(config["exacerbation_severity"]),
            antibiotic_exposure=AntibioticExposure(config["antibiotic_exposure"], abx_mid_trends, nothing),
            family_history=FamilyHistory(config["family_history"]),
            utility=Utility(config["utility"]),
            cost=AsthmaCost(config["cost"]),
            census_table=CensusTable(config["census_table"]),
            initial_distribution=nothing,
            outcome_matrix=(;)
        )

        return simulation
    else
        error("Province not supported")
    end
end
