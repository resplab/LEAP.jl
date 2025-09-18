using CSV, DataFrames, Logging
using Roots
include("utils.jl")

DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
STARTING_YEAR = 1996
FINAL_YEAR = 2068

DESIRED_LIFE_EXPECTANCIES = DataFrame(
    :province => Array{AbstractString, 1}(["CA", "CA", "BC", "BC"]),
    :sex => Array{AbstractString, 1}(["M", "F", "M", "F"]),
    :life_expectancy => Array{Float64, 1}([87.0, 90.1, 84.6, 88.0])
)

CALIBRATION_YEARS = Dict(
    "CA" => 2068,
    "BC" => 2043,
    "AB" => missing,
    "SK" => missing,
    "MB" => missing,
    "ON" => missing,
    "QC" => missing,
    "NL" => missing,
    "NS" => missing,
    "NB" => missing,
    "PE" => missing,
    "YT" => missing,
    "NT" => missing,
    "NU" => missing
)


function calculate_life_expectancy(life_table::DataFrame)::Float64
    """
        calculate_life_expectancy(life_table)

    Determine the life expectancy for a person born in a given year.

    The life expectancy can be calculated from the death probability using the formulae
    delineated here: https://www.ssa.gov/oact/HistEst/CohLifeTables/LifeTableDefinitions.pdf.

    # Arguments
    - `life_table::DataFrame`: A dataframe containing the probability of death for a single year,
        province and sex, for each age. Columns:
            - `age`: the integer age.
            - `sex`: one of "M" or "F".
            - `year`: the integer calendar year.
            - `province`: a string indicating the province abbreviation, e.g. "BC".
                For all of Canada, set province to "CA".
            - `prob_death`: the probability of death for a given age, province, sex, and year.

    # Returns
    - `Float64`: The life expectancy for a person born in the given year, in a given province,
        for a given sex.
    """
    n_alive_by_age_0 = 100000
    n_alive_by_age = Array{Float64, 1}()
    for age in life_table.age
        age_index = age + 1
        if age == 0
            push!(n_alive_by_age, n_alive_by_age_0)
        else
            push!(
                n_alive_by_age,
                n_alive_by_age[age_index - 1] * (1 - life_table[age_index - 1, :prob_death])
            )
        end
    end
    life_table[!, :n_alive_by_age] = n_alive_by_age

    life_table = transform(
        life_table,
        [:prob_death, :n_alive_by_age] =>
            ByRow((x, y) -> y - 0.5 * x * y) => :n_person_years_interval
    )

    life_table[1, :n_person_years_interval] = (
        life_table[2, :n_person_years_interval] +
        0.1 * life_table[1, :prob_death] * n_alive_by_age_0
    )

    life_table[111, :n_person_years_interval] = life_table[111, :n_alive_by_age] * 1.4

    life_table = select(
        life_table,
        :,
        :n_person_years_interval => ((x) -> sum(x)) => :n_person_years_after_age
    )

    life_table = transform(
        life_table,
        [:n_person_years_after_age, :n_alive_by_age] =>
            ByRow((x, y) -> x / y) => :n_years_left_at_age
    )

    life_expectancy = life_table[1, :n_years_left_at_age]
    return life_expectancy
end


function get_prob_death_projected(
    prob_death::Float64, year_index::Integer, beta_year::Float64
)::Float64
    """
        get_prob_death_projected(prob_death, year_index, beta_year)

    Given the probability of death for a past year, calculate the probability of death in a
    future year.

    # Arguments
    - `prob_death::Float64`: The probability of death for the initial year
        (determined by past data).
    - `year_index::Integer`: The number of years between the current year and the initial year.
        For example, if our initial year is 2020, and we want to compute the probability of death
        in 2028, the `year_index` would be 8.
    - `beta_year::Float64`: The beta parameter for the given year.

    # Returns
    - `Float64`: The projected probability of death for the current year.
    """
    prob_death = min(prob_death, 0.9999999999)
    odds = (prob_death / (1 - prob_death)) * exp(year_index * beta_year)
    prob_death_projected = max(min(odds / (1 + odds), 1), 0)
    return prob_death_projected
end


function get_projected_life_table_single_year(
    beta_year::Float64, life_table::DataFrame, starting_year::Integer,
    year_index::Integer, sex::AbstractString, province::AbstractString
)::DataFrame
    """
        get_projected_life_table_single_year(
            beta_year, life_table, starting_year, year_index, sex, province
        )

    Given the probability of death for a past year, calculate the probability of death in a
    future year.

    # Arguments
    - `beta_year::Float64`: The beta parameter for the given year.
    - `life_table::DataFrame`: A dataframe containing the projected probability of death
        for the starting year, for a given sex and province. Columns:
            - `age`: the integer age.
            - `sex`: one of "M" or "F".
            - `year`: the starting calendar year.
            - `province`: a string indicating the province abbreviation, e.g. "BC".
                For all of Canada, set province to "CA".
            - `prob_death`: the probability of death for a given age, province, sex, and year.
    - `starting_year::Integer`: The calendar year when the projections begin.
    - `year_index::Integer`: The number of years between the current year and the starting year.
        For example, if our initial year is 2020, and we want to compute the probability of death
        in 2028, the `year_index` would be 9.
    - `sex::AbstractString`: one of "M" or "F".
    - `province::AbstractString`: a string indicating the province abbreviation, e.g. "BC".
        For all of Canada, set province to "CA".

    # Returns
    - `DataFrame`: A dataframe containing the projected probability of death for the given year,
        sex, and province.
    """

    df = filter(([:sex] => (x) -> x == sex), life_table)
    df = transform(
        df,
        [:prob_death] =>
            ByRow((x) -> get_prob_death_projected(x, year_index, beta_year))
            => :prob_death_proj
    )
    df = transform(
        df,
        [:year] =>
            ByRow((x) -> starting_year + year_index - 1)
            => :year
    )
    df = transform(
        df,
        [:se, :prob_death, :prob_death_proj] =>
        ByRow((x, y, z) -> (x * z) / y)
        => :se
    )

    df = select(df, Not([:prob_death]))
    df = rename(df, [:prob_death_proj => :prob_death])

    return df
end


function beta_year_optimizer(beta_year::Float64, p::AbstractDict)::Float64
    """
        beta_year_optimizer(beta_year, p)

    Calculate the difference between the projected life expectancy and the desired life expectancy.

    This function is passed to the `Roots.find_zero` function. We want to find `beta_year` such
    that the projected life expectancy is as close as possible to the desired life expectancy.
    Unfortunately, the `Roots.find_zero` function does not have the most user-friendly options
    to pass arguments to the function, so we have to pass all additional arguments as a dictionary
    called `p`.

    # Arguments
    - `beta_year::Float64`: The beta parameter for the given year.
    - `p::AbstractDict`: A dictionary with the following keys:
        - `year_index::Integer`: The number of years between the current year and the starting year.
            For example, if our initial year is 2020, and we want to compute the probability of
            death in 2028, the `year_index` would be 9.
        - `sex`: one of "M" or "F".
        - `province`: a string indicating the province abbreviation, e.g. "BC".
            For all of Canada, set province to "CA".
        - `life_table::DataFrame`: A dataframe containing the projected probability of death
            for the calibration year, for a given sex and province. Columns:
            - `age`: the integer age.
            - `sex`: one of "M" or "F".
            - `year`: the calibration calendar year.
            - `province`: a string indicating the province abbreviation, e.g. "BC".
                For all of Canada, set province to "CA".
            - `prob_death`: the probability of death for a given age, province, sex, and year.

    # Returns
    - `Float64`: The difference between the projected life expectancy of the calibration year
        and the desired life expectancy.
    """

    year_index = p[:year_index]
    sex = p[:sex]
    province = p[:province]
    life_table = p[:life_table]
    starting_year = p[:starting_year]

    projected_life_table = get_projected_life_table_single_year(
        beta_year, life_table, starting_year, year_index, sex, province
    )

    life_expectancy = calculate_life_expectancy(projected_life_table)
    desired_life_expectancy = filter(
        ([:sex, :province] => (x, y) ->
            x == sex && y == province
        ),
        DESIRED_LIFE_EXPECTANCIES
    )[1, :life_expectancy]

    return life_expectancy - desired_life_expectancy
end


function load_past_death_data()
    @info "Loading mortality data from CSV file..."
    df = CSV.read(
        joinpath(DATA_PATH, "13100837.csv"),
        DataFrame
    )

    # remove spaces from column names and make uppercase
    column_names = Array{Symbol, 1}()
    for column in names(df)
        push!(column_names, Symbol(uppercase(replace(column, " " => "_"))))
    end
    df = rename(df, column_names)

    # rename the columns
    df = rename(
        df,
        [
            :REF_DATE => :year, :GEO => :province, :SEX => :sex, :AGE_GROUP => :age
        ]
    )

    # select the required columns
    df = select(
        filter(([:year] => (x) -> x >= STARTING_YEAR), df),
        [:year, :province, :sex, :age, :ELEMENT, :VALUE]
    )

    # format the age group into an integer age
    df = transform(df, :age => ByRow(x -> format_age_group(x, "110 years and over")) => :age)

    # convert province names to 2-letter province IDs
    df = transform(df, :province => ByRow(x -> get_province_id(x)) => :province)

    # filter only "CA" and "BC"
    df = filter(([:province] => (x) -> x in ["CA", "BC"]), df)

    # convert :sex to 1-letter ID ("F", "M", "B")
    df = transform(df, :sex => ByRow(x -> get_sex_id(x)) => :sex)

    # remove :sex category "Both"
    df = filter(([:sex] => (x) -> x != "B"), df)

    # select only the "qx" elements, which relate to the probability of death and the SE
    df = filter(([:ELEMENT] => (x) -> occursin("qx", x)), df)

    # create a df with the probability of death
    df_prob = filter(([:ELEMENT] => (x) -> occursin("Death probability between age", x)), df)
    df_prob = select(df_prob, Not(:ELEMENT))
    df_prob = rename(df_prob, [:VALUE => :prob_death])

    # create a df with the standard error of the probability of death
    df_se = filter(([:ELEMENT] => (x) -> occursin("Margin of error", x)), df)
    df_se = select(df_se, Not(:ELEMENT))
    df_se = rename(df_se, [:VALUE => :se])

    # join the two tables
    df = leftjoin(df_prob, df_se, on=[:year, :province, :sex, :age])
    return df
end


function load_projected_death_data(past_life_table)

    projected_life_table = DataFrame(
        :year => Integer[],
        :province => AbstractString[],
        :age => Integer[],
        :sex => AbstractString[],
        :prob_death => Float64[],
        :se => Float64[]
    )
    for province in Set(past_life_table[!, :province])
        calibration_year = CALIBRATION_YEARS[province]
        life_table = filter(([:province] => (x) -> x == province), past_life_table)
        starting_year = findmax(life_table[!, :year])[1] + 1
        life_table = filter(([:year] => (x) -> x == starting_year - 1), life_table)
        n_years = FINAL_YEAR - starting_year + 1

        beta_year_female = find_zero(
            beta_year_optimizer,
            (-0.03, -0.01),
            Roots.Brent(),
            p = Dict(
                :life_table => life_table,
                :sex => "F",
                :province => province,
                :year_index => calibration_year - starting_year + 1,
                :starting_year => starting_year
            ),
            atol=0.00001
        )
        beta_year_male = find_zero(
            beta_year_optimizer,
            (-0.03, -0.01),
            Roots.Brent(),
            p=Dict(
                :life_table => life_table,
                :sex => "M",
                :province => province,
                :year_index => calibration_year - starting_year + 1,
                :starting_year => starting_year
            ),
            atol=0.00001
        )
        @info beta_year_male
        projected_life_table_province = DataFrame(
            :year => Integer[],
            :province => AbstractString[],
            :age => Integer[],
            :sex => AbstractString[],
            :prob_death => Float64[],
            :se => Float64[]
        )
        for year_index in 1:n_years
            df_female = get_projected_life_table_single_year(
                beta_year_female, life_table, starting_year, year_index, "F", province
            )
            df_male = get_projected_life_table_single_year(
                beta_year_male, life_table, starting_year, year_index, "M", province
            )
            projected_life_table_single_year = vcat(df_female, df_male)
            append!(projected_life_table_province, projected_life_table_single_year)
        end

        append!(projected_life_table, projected_life_table_province)
    end
    return projected_life_table
end


function process_death_data()
    past_life_table = load_past_death_data()
    projected_life_table = load_projected_death_data(past_life_table)
    life_table = vcat(past_life_table, projected_life_table)
    # save the data
    file_path = joinpath(PROCESSED_DATA_PATH, "master_life_table.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, life_table)
end

process_death_data()
