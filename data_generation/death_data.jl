using CSV, DataFrames, Logging
include("utils.jl")

DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
STARTING_YEAR = 1996


function load_death_data()
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

    # save the data
    file_path = joinpath(PROCESSED_DATA_PATH, "life_table.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, df)
end
