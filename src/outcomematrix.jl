"""
    OutcomeMatrix

TODO.

# Fields
- `immigration::Union{Array{Int}, Nothing}`: TODO.
"""
@kwdef mutable struct OutcomeMatrix  <: OutcomeMatrixModule
    asthma_incidence_contingency_table::Union{GroupedDataFrame{DataFrame}, Matrix{Int}, Nothing}
    asthma_prevalence_contingency_table::Union{GroupedDataFrame{DataFrame}, Matrix{Int}, Nothing}
    control::Union{Array{Real}, Nothing}
    cost::Union{Array{Real}, Nothing}
    utility::Union{Array{Real}, Nothing}
    exacerbation_by_severity::Union{Array{Real}, Nothing}
    asthma_incidence_family_history::Union{Array{Int}, Nothing}
    asthma_prevalence_family_history::Union{Array{Int}, Nothing}
    asthma_status_family_history::Union{Array{Int}, Nothing}
    family_history_prevalence::Union{Array{Int}, Nothing}
    asthma_incidence_antibiotic_exposure::Union{Array{Int}, Nothing}
    asthma_prevalence_antibiotic_exposure::Union{Array{Int}, Nothing}
    asthma_status_antibiotic_exposure::Union{Array{Int}, Nothing}
    asthma_status::Union{Array{Int}, Nothing}
    antibiotic_exposure::Union{Array{Int}, Nothing}
    asthma_incidence::Union{Array{Int}, Nothing}
    asthma_prevalence::Union{Array{Int}, Nothing}
    death::Union{Array{Int}, Nothing}
    alive::Union{Array{Int}, Nothing}
    family_history::Union{Array{Int}, Nothing}
    exacerbation::Union{Array{Int}, Nothing}
    exacerbation_hospital::Union{Array{Int}, Nothing}
    emigration::Union{Array{Int}, Nothing}
    immigration::Union{Array{Int}, Nothing}
end


function set_outcome_matrix_asthma_prevalence_contingency_table!(
    outcome_matrix::OutcomeMatrix,
    cal_years::Union{UnitRange{Int}, Nothing}=nothing,
    min_year::Union{Integer, Nothing}=nothing,
    max_year::Union{Integer, Nothing}=nothing,
    until_all_die::Union{Bool, Nothing}=nothing,
    max_age::Union{Integer, Nothing}=nothing,
    asthma_prevalence_contingency_table::Union{GroupedDataFrame{DataFrame}, Nothing}=nothing
)
    outcome_matrix.asthma_prevalence_contingency_table = create_contingency_table(
        cal_years, min_year, max_year, until_all_die, max_age
    )
    return outcome_matrix
end

function set_outcome_matrix_asthma_incidence_contingency_table!(
    outcome_matrix::OutcomeMatrix,
    cal_years::Union{UnitRange{Int}, Nothing}=nothing,
    min_year::Union{Integer, Nothing}=nothing,
    max_year::Union{Integer, Nothing}=nothing,
    until_all_die::Union{Bool, Nothing}=nothing,
    max_age::Union{Integer, Nothing}=nothing,
    asthma_incidence_contingency_table::Union{GroupedDataFrame{DataFrame}, Nothing}=nothing
)
    outcome_matrix.asthma_incidence_contingency_table = create_contingency_table(
        cal_years, min_year, max_year, until_all_die, max_age
    )
    return outcome_matrix
end

function create_contingency_table(
    cal_years::Union{UnitRange{Int}, Nothing}=nothing,
    min_year::Union{Integer, Nothing}=nothing,
    max_year::Union{Integer, Nothing}=nothing,
    until_all_die::Union{Bool, Nothing}=nothing,
    max_age::Union{Integer, Nothing}=nothing,
    contingency_table::Union{GroupedDataFrame{DataFrame}, Nothing}=nothing
)
    if isnothing(contingency_table)
        tmp_df = DataFrame(year=Int64[], sex=Int64[], age=Int64[], fam_history = Int64[],
            abx_exposure = Int64[], n_asthma= Int64[], n_no_asthma = Int64[])
        foreach(x -> push!(tmp_df, x), Iterators.product(
            min_year:1:max_year, 0:1:1, 0:1:max_age+1, 0:1:1, 0:1:3, 0, 0)
        )
        tmp_df = groupby(tmp_df, [:year, :sex, :fam_history, :abx_exposure])
        return tmp_df
    else
        return contingency_table
    end
end


function setfield_outcome_matrix!(
    outcome_matrix::OutcomeMatrix,
    field::String,
    type::DataType,
    dimensions::Tuple{Vararg{Int, N}},
    df::Union{GroupedDataFrame{DataFrame}, Nothing}=nothing
) where N
    if isnothing(df)
        setproperty!(outcome_matrix, Symbol(field), zeros(type, dimensions))
    else
        setproperty!(outcome_matrix, Symbol(field), df)
    end
    return outcome_matrix
end



function create_outcome_matrix(; until_all_die::Bool, cal_years::UnitRange{Int},
    min_year::Integer, max_year::Integer, max_age::Integer)

    outcome_matrix = OutcomeMatrix(
        control=nothing, cost=nothing, utility=nothing, exacerbation_by_severity=nothing,
        asthma_incidence_family_history=nothing,
        asthma_prevalence_family_history=nothing,
        asthma_status_family_history=nothing,
        family_history_prevalence=nothing,
        asthma_incidence_antibiotic_exposure=nothing,
        asthma_prevalence_antibiotic_exposure=nothing,
        asthma_status_antibiotic_exposure=nothing,
        asthma_status=nothing,
        antibiotic_exposure=nothing,
        asthma_incidence=nothing,
        asthma_prevalence=nothing,
        family_history=nothing,
        death=nothing,
        alive=nothing,
        exacerbation=nothing,
        exacerbation_hospital=nothing,
        emigration=nothing,
        immigration=nothing,
        asthma_incidence_contingency_table=nothing,
        asthma_prevalence_contingency_table=nothing
    )

    set_outcome_matrix_asthma_incidence_contingency_table!(
        outcome_matrix, cal_years, min_year, max_year, until_all_die, max_age
    )
    set_outcome_matrix_asthma_prevalence_contingency_table!(
        outcome_matrix, cal_years, min_year, max_year, until_all_die, max_age
    )

    type = Real
    dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2, 3)
    setfield_outcome_matrix!(outcome_matrix, "control", type, dimensions)

    type = Real
    dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2, 4)
    setfield_outcome_matrix!(outcome_matrix, "exacerbation_by_severity", type, dimensions)

    type = Int
    dimensions = (2, length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
    setfield_outcome_matrix!(outcome_matrix, "asthma_incidence_family_history", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_prevalence_family_history", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_status_family_history", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "family_history_prevalence", type, dimensions)

    type = Int
    dimensions = (4, length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
    setfield_outcome_matrix!(outcome_matrix, "asthma_incidence_antibiotic_exposure", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_prevalence_antibiotic_exposure", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_status_antibiotic_exposure", type, dimensions)

    type = Real
    dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
    setfield_outcome_matrix!(outcome_matrix, "cost", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "utility", type, dimensions)

    type = Int
    dimensions = (length(cal_years) + (until_all_die ? max_age : 0), max_age + 1, 2)
    setfield_outcome_matrix!(outcome_matrix, "asthma_status", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "antibiotic_exposure", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_incidence", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "asthma_prevalence", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "death", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "alive", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "exacerbation", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "exacerbation_hospital", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "emigration", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "immigration", type, dimensions)
    setfield_outcome_matrix!(outcome_matrix, "family_history", type, dimensions)

    return outcome_matrix
end


function add_control_to_outcome_matrix!(outcome_matrix::OutcomeMatrix, age::Integer, sex::Bool,
    cal_year_index::Integer, control_levels::AbstractDict)
    outcome_matrix.control[cal_year_index, age + 1, sex + 1, :] += control_levels[:as_array]
end

function add_exacerbation_by_severity_to_outcome_matrix!(outcome_matrix::OutcomeMatrix, age::Integer,
    sex::Bool, cal_year_index::Integer,
    exac_sev_hist_current_year::Vector{Integer})
    outcome_matrix.exacerbation_by_severity[cal_year_index, age + 1, sex + 1, :] .+= exac_sev_hist_current_year
end

function update_asthma_in_contingency_table!(outcome_matrix::OutcomeMatrix,
    age::Integer, sex::Bool, cal_year::Integer, has_family_hist::Bool, num_antibiotic_use::Integer,
    has_asthma::Bool, inc_or_prev::String
)
    if has_asthma
        column = "n_asthma"
    else
        column = "n_no_asthma"
    end

    if inc_or_prev == "incidence"
        outcome_matrix.asthma_incidence_contingency_table[(
            cal_year, Int(sex),
            Int(has_family_hist),
            min(num_antibiotic_use, 3)
            )][age + 1, column] += 1
    elseif inc_or_prev == "prevalence"
        outcome_matrix.asthma_prevalence_contingency_table[(
            cal_year, Int(sex),
            Int(has_family_hist),
            min(num_antibiotic_use, 3)
            )][age + 1, column] += 1
    else
        throw(ArgumentError(
            "inc_or_prev must be either 'prevalence' or 'incidence', received $inc_or_prev."
        ))
    end
end



function increment_field_in_outcome_matrix!(outcome_matrix::OutcomeMatrix, field::String,
    age::Integer, sex::Bool, cal_year_index::Integer, increment::Float64=1.0)

    field_object = getfield(outcome_matrix, Symbol(field))
    field_object[cal_year_index, age + 1, sex + 1] += increment
    setproperty!(outcome_matrix, Symbol(field), field_object)
end

function increment_field_in_outcome_matrix!(outcome_matrix::OutcomeMatrix, field::String,
    age::Integer, sex::Bool, cal_year_index::Integer, increment::Integer=1)

    field_object = getfield(outcome_matrix, Symbol(field))
    field_object[cal_year_index, age + 1, sex + 1] += increment
    setproperty!(outcome_matrix, Symbol(field), field_object)
end
