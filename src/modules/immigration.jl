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


"""
    get_num_new_immigrants(immigration, num_new_born, cal_year_index)

Get the number of new immigrant to Canada in a given year.

# Arguments
- `immigration::Immigration`: an Immigration object, see  [`Immigration`](@ref).
- `num_new_born::Integer`: number of births in the given year of the simulation.
- `cal_year_index::Integer`: An integer representing the year of the simulation. For example, if the
    simulation starts in 2023, then the `cal_year_index` for 2023 is 1, for 2024 is 2, etc.

# Returns
- `Integer`: the number of new immigrants to Canada in a given year.
"""
function get_num_new_immigrants(immigration::Immigration, num_new_born::Integer, cal_year_index::Integer)
    num_new_immigrants = ceil(
        Int,
        num_new_born * sum(immigration.table[cal_year_index].n_prop_birth)
    )
    return num_new_immigrants
end
