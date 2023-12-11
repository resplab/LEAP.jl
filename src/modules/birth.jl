struct Birth <: Birth_Module
    estimate
    initial_population
end

function process_birth(cal_year::Integer, cal_year_index::Integer, birth::Birth)
    agent = Agent(
        sex=rand(Bernoulli(birth.estimate.prop_male[cal_year_index])),
        age=0, cal_year=cal_year, cal_year_index=cal_year_index, alive=true, num_antibiotic_use=0,
        has_asthma=false, asthma_age=nothing, severity=nothing, control=nothing,
        exac_hist=ExacerbationHist(0, 0), exac_sev_hist=[zeros(4),zeros(4)], total_hosp=0,
        family_hist=false, asthma_status=false
    )
    return agent
end


function process_birth(cal_year::Integer, cal_year_index::Integer, birth::Birth,
    antibiotic_exposure::AntibioticExposure, family_hist::FamilyHistory)
    agent = Agent(
        sex=rand(Bernoulli(birth.estimate.prop_male[cal_year_index])),
        age=0, cal_year=cal_year, cal_year_index=cal_year_index, alive=true, num_antibiotic_use=0,
        has_asthma=false, asthma_age=nothing, severity=nothing, control=nothing,
        exac_hist=ExacerbationHist(0, 0), exac_sev_hist=[zeros(4),zeros(4)],
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


function process_birth(cal_year::Integer, cal_year_index::Integer, birth::Birth, sex::Bool,
    age::Integer, antibiotic_exposure::AntibioticExposure, family_hist::FamilyHistory)
    agent = Agent(
        sex=sex, age=age, cal_year=cal_year, cal_year_index=cal_year_index, alive=true,
        num_antibiotic_use=0, has_asthma=false, asthma_age=nothing, severity=nothing,
        control=nothing, exac_hist=ExacerbationHist(0, 0),
        exac_sev_hist=[zeros(4),zeros(4)], total_hosp=0,
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
