"""
    Reassessment

A struct containing information about asthma diagnosis reassessment.

# Fields
- `table::GroupedDataFrame{DataFrame}`: A grouped data frame grouped by year.
    Each data frame contains the following columns:
        `year::Integer`: year.
        `age::Integer`: age of person.
        `F::Float64`: the probability that a female agent still has asthma.
        `M::Float64`: the probability that a male agent still has asthma.
        `province::String`: a string indicating the province abbreviation, e.g. "BC".
            For all of Canada, set province to "CA".
    See `master_asthma_reassessment.csv`.
"""
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


"""
    load_reassessment_table(starting_year, province)

Load the asthma diagnosis reassessment table.

# Arguments
- `starting_year::Integer`: the year to start the data at.
- `province::String`: a string indicating the province abbreviation, e.g. "BC". For all of Canada,
    set province to "CA".

# Returns
- `GroupedDataFrame{DataFrame}`: A grouped data frame grouped by year.
    Each data frame contains the following columns:
        `year::Integer`: year.
        `age::Integer`: age of person.
        `F::Float64`: TODO.
        `M::Float64`: TODO.
        `province::String`: a string indicating the province abbreviation, e.g. "BC".
            For all of Canada, set province to "CA".
"""
function load_reassessment_table(starting_year::Integer, province::String)
    master_reassessment = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_asthma_reassessment.csv"),
        DataFrame
    )
    reassessment_table = groupby(
        filter(
            [:year, :province] => (x, y) ->
            x >= starting_year &&
            y == province,
            master_reassessment
        ),
        :year
    )
    return reassessment_table
end


"""
    agent_has_asthma(agent, reassessment)

If an agent has been diagnosed with asthma, determine whether the agent still has asthma.

Asthma is not curable, but it can be mistaken for other respiratory diseases, and so a person could
be diagnosed with asthma and later find this was a misdiagnosis. Additionally, asthma can go into
a period of dormancy, during which time the person does not have any asthma symptoms.

# Arguments

- `agent::Agent`: Agent module, see [`Agent`](@ref).
- `reassessment::Reassessment`: Reassessment module, see [`Reassessment`](@ref).
"""
function agent_has_asthma(agent::Agent, reassessment::Reassessment)
    max_year = length(reassessment.table)
    if agent.age < 4
        agent.has_asthma
    else
        rand(Bernoulli(
            reassessment.table[min(agent.cal_year_index, max_year)][agent.age - 3, agent.sex + 3]
        ))
    end
end
