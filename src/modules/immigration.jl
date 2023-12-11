# include("abstractModule.jl")
# include("../utils.jl")

# estimated num of newborn + proprotion of male
# birth_project = CSV.read("../processed_data/brith_projection.csv",DataFrame)

struct Immigration <: Immigration_Module
    sex_ratio
    estimate
    age_distribution
    overall_rate
    table
end

function process_immigration(sex::Bool, age::Integer, cal_year::Integer, cal_year_index::Integer)
    agent = Agent(
        sex=sex, age=age, cal_year=cal_year, cal_year_index=cal_year_index, alive=true,
        num_antibiotic_use=0, has_asthma=false, asthma_age=nothing, severity=nothing,
        control=nothing, exac_hist=ExacerbationHist(0, 0),
        exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)), total_hosp=0,
        family_hist=false, asthma_status=false
    )
    return agent
end


function process_immigration(sex::Bool, age::Integer, cal_year::Integer, cal_year_index::Integer,
    antibiotic_exposure::AntibioticExposure, family_hist::FamilyHistory)
    agent = Agent(
        sex=sex, age=age, cal_year=cal_year, cal_year_index=cal_year_index, alive=true,
        num_antibiotic_use=0, has_asthma=false, asthma_age=nothing, severity=nothing,
        control=nothing, exac_hist=ExacerbationHist(0, 0),
        exac_sev_hist=ExacerbationSeverityHist(zeros(4),zeros(4)), total_hosp=0,
        family_hist=false, asthma_status=false
    )
    @set! agent.num_antibiotic_use = process_initial(
        agent, antibiotic_exposure, cal_year - age
    )
    @set! agent.family_hist = process_initial(agent, family_hist)
    return agent
end
