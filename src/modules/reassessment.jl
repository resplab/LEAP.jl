struct Reassessment <: ReassessmentModule
    table::GroupedDataFrame{DataFrame}
    function Reassessment(starting_year::Integer, province::String)
        table = load_reassessment_table(starting_year, province)
        new(table)
    end
    function Reassessment(table::GroupedDataFrame{DataFrame})
        new(table)
    end
end


function load_reassessment_table(starting_year::Integer, province::String)
    master_reassessment = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_assessment.csv"),
        DataFrame
    )
    reassessment_table = groupby(
        filter(
            [:year, :province] => (x, y) -> x >= starting_year && y == province,
            master_reassessment
        ),
        :year
    )
    return reassessment_table
end


function agent_has_asthma(agent::Agent, reassessment::Reassessment)
    max_year = length(reassessment.table)
    if agent.age < 4
        agent.has_asthma
    else
        rand(Bernoulli(reassessment.table[min(agent.cal_year_index,max_year)][agent.age-3, agent.sex+3]))
    end
end
