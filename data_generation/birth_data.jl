using CSV, DataFrames, Logging
include("utils.jl")

DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
STARTING_YEAR = 1999


function get_projection_scenario_id(projection_scenario)
    projection_scenario_id = replace(projection_scenario, "Projection scenario " => "")
    return projection_scenario_id[1:2]
end

function filter_age_group(age_group::AbstractString)
    if occursin("100", age_group)
        return true
    else
        return (
            !occursin("to", age_group) && !occursin("over", age_group)
            && !occursin("All", age_group) && !occursin("Median", age_group)
            && !occursin("Average", age_group)
        )
    end
end


function load_past_population_data()
    @info "Loading past population data from CSV file..."
    df = CSV.read(
        joinpath(DATA_PATH, "17100005.csv"),
        DataFrame
    )

    df = select(
        filter(
            ([:REF_DATE, :AGE_GROUP] => (x, y) ->
                x >= STARTING_YEAR &&
                y == "0 years"
            ),
            df
        ),
        [:REF_DATE, :GEO, :SEX, :VALUE]
    )
    df = rename(df, [:REF_DATE => :year, :GEO => :province, :SEX => :sex, :VALUE => :N])

    # convert province names to 2-letter province IDs
    df = transform(df, :province => ByRow(x -> get_province_id(x)) => :province)

    # convert :sex to 1-letter ID ("F", "M", "B")
    df = transform(df, :sex => ByRow(x -> get_sex_id(x)) => :sex)

    # convert :N to integer
    df = transform(df, :N => ByRow(x -> Int(x)) => :N)

    # get the proportion male / female
    df = groupby(df, [:year, :province])
    df = select(df, :, :N => ((x) -> findmax(x)[1]) => :max_N, ungroup=false)
    df = transform(df, [:N, :max_N] => ByRow((x, y) -> x / y) => :prop, ungroup=true)

    # keep only male entries
    df = filter(([:sex] => (x) -> x == "M"), df)

    # drop :N and :sex columns
    df = select(df, Not(:N, :sex))
    df = rename(df, [:max_N => :N, :prop => :prop_male])

    # add :projection_scenario column, all values = "past"
    df[:, "projection_scenario"] = repeat(["past"], size(df)[1])

    return df
end


function load_projected_population_data(min_year::Integer)
    @info "Loading projected population data from CSV file..."
    df = CSV.read(
        joinpath(DATA_PATH, "17100057.csv"),
        DataFrame
    )

    # remove spaces from column names and make uppercase
    column_names = Array{Symbol, 1}()
    for column in names(df)
        push!(column_names, Symbol(uppercase(replace(column, " " => "_"))))
    end
    df = rename(df, column_names)

    df = select(
        filter(
            ([:REF_DATE, :AGE_GROUP] => (x, y) -> x >= min_year && y == "Under 1 year"),
            df
        ),
        [:REF_DATE, :GEO, :PROJECTION_SCENARIO, :SEX, :AGE_GROUP, :VALUE]
    )

    df = rename(
        df,
        [
            :REF_DATE => :year, :GEO => :province, :SEX => :sex, :AGE_GROUP => :age,
            :VALUE => :N, :PROJECTION_SCENARIO => :projection_scenario
        ]
    )

    # convert the long form of the projection scenario to the 2-letter ID
    transform(
        df,
        :projection_scenario => ByRow(x -> get_projection_scenario_id(x)) => :projection_scenario
    )

    # convert province names to 2-letter province IDs
    df = transform(df, :province => ByRow(x -> get_province_id(x)) => :province)

    # convert :sex to 1-letter ID ("F", "M", "B")
    df = transform(df, :sex => ByRow(x -> get_sex_id(x)) => :sex)

    # format the age group string
    df = transform(df, :age => ByRow(x -> 0) => :age)

    # remove rows which are missing values of N
    df = dropmissing(df, :N)

    # multiply the :N column by 1000 and convert to integer
    df = transform(df, :N => ByRow(x -> Int(round(x * 1000, digits=0))) => :N)

    # get the proportion male / female
    df = groupby(df, [:year, :province, :projection_scenario])
    df = select(df, :, :N => ((x) -> findmax(x)[1]) => :max_N, ungroup=false)
    df = transform(df, [:N, :max_N] => ByRow((x, y) -> x / y) => :prop, ungroup=true)

    # keep only male entries
    df = filter(([:sex] => (x) -> x == "M"), df)

    # drop :N and :sex columns
    df = select(df, Not(:N, :sex, :age))
    df = rename(df, [:max_N => :N, :prop => :prop_male])
    return df
end


function load_past_initial_population_data()
    @info "Loading past population data from CSV file..."
    df = CSV.read(
        joinpath(DATA_PATH, "17100005.csv"),
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
            :REF_DATE => :year, :GEO => :province, :SEX => :sex, :AGE_GROUP => :age, :VALUE => :N
        ]
    )

    # select the required columns
    df = select(
        filter(([:year] => (x) -> x >= STARTING_YEAR + 1), df),
        [:year, :province, :sex, :age, :N]
    )

    df = filter(([:age] => (x) -> filter_age_group(x)), df)
    df = transform(df, :age => ByRow(x -> format_age_group(x)) => :age)

    # convert province names to 2-letter province IDs
    df = transform(df, :province => ByRow(x -> get_province_id(x)) => :province)

    # convert :sex to 1-letter ID ("F", "M", "B")
    df = transform(df, :sex => ByRow(x -> get_sex_id(x)) => :sex)

    # remove :sex category "Both"
    df = filter(([:sex] => (x) -> x != "B"), df)

    # find the missing values of :N
    missing_df = filter(([:N] => (x) -> ismissing(x)), df)
    missing_df = select(missing_df, Not(:N))

    # create a df to replace missing values with those of the next year and age
    replacement_df = filter(
        ([:year, :age] => (x, y) ->
            x in Set(missing_df.year .+ 1) &&
            y in Set(missing_df.age .+ 1)
        ),
        df
    )
    replacement_df = transform(replacement_df, :age => ByRow(x -> x - 1) => :age)
    replacement_df = select(replacement_df, Not(:year))
    replacement_df = rename(replacement_df, [:N => :N_replace])

    # merge the two dfs
    replacement_df = leftjoin(missing_df, replacement_df, on=[:sex, :age, :province])

    # replace the missing values in the original df
    df = leftjoin(df, replacement_df, on=[:sex, :age, :province, :year])
    df = transform(df, [:N, :N_replace]  => ByRow((x, y) -> ismissing(x) ? y : x) => :N)
    df = select(df, Not(:N_replace))

    # remove rows which are still missing values of N
    df = dropmissing(df, :N)

    # convert :N to integer
    df = transform(df, :N => ByRow(x -> Int(round(x, digits=0))) => :N)

    # get the total population for a given year, province, and age
    df = groupby(df, [:year, :age, :province])
    df = select(df, :, :N => ((x) -> sum(x)) => :N_age, ungroup=false)
    df = transform(df, [:N, :N_age] => ByRow((x, y) -> x / y) => :prop_male, ungroup=true)

    # get the total number of births for a given year and province
    df_birth = filter(([:age] => (x) -> x == 0), df)
    df_birth = transform(df_birth, [:N_age] => ByRow((x) -> x) => :N_birth)
    df_birth = select(df_birth, Not([:age, :N, :N_age, :prop_male]))

    # add the births column to the main df
    df = leftjoin(df, df_birth, on=[:province, :sex, :year])
    df = transform(df, [:N_age, :N_birth] => ByRow((x, y) -> x / y) => :prop)

    # keep only male entries
    df = filter(([:sex] => (x) -> x == "M"), df)
    df = select(df, Not([:sex]))

    # add :projection_scenario column, all values = "past"
    df[:, "projection_scenario"] = repeat(["past"], size(df)[1])

    return df
end


function load_projected_initial_population_data(min_year::Integer)
    @info "Loading projected population data from CSV file..."
    df = CSV.read(
        joinpath(DATA_PATH, "17100057.csv"),
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
            :REF_DATE => :year, :GEO => :province, :SEX => :sex, :AGE_GROUP => :age,
            :VALUE => :N, :PROJECTION_SCENARIO => :projection_scenario

        ]
    )

    # select the required columns
    df = select(
        filter(([:year] => (x) -> x >= min_year), df),
        [:year, :province, :sex, :age, :N, :projection_scenario]
    )

    # convert the long form of the projection scenario to the 2-letter ID
    df = transform(
        df,
        :projection_scenario => ByRow(x -> get_projection_scenario_id(x)) => :projection_scenario
    )

    df = filter(([:age] => (x) -> filter_age_group(x)), df)
    df = transform(df, :age => ByRow(x -> format_age_group(x)) => :age)

    # convert province names to 2-letter province IDs
    df = transform(df, :province => ByRow(x -> get_province_id(x)) => :province)

    # convert :sex to 1-letter ID ("F", "M", "B")
    df = transform(df, :sex => ByRow(x -> get_sex_id(x)) => :sex)

    # remove :sex category "Both"
    df = filter(([:sex] => (x) -> x != "B"), df)

    # remove rows which are missing values of N
    df = dropmissing(df, :N)

    # multiply the :N column by 1000 and convert to integer
    df = transform(df, :N => ByRow(x -> Int(round(x * 1000, digits=0))) => :N)

    # get the total population for a given year, province, age, and projection scenario
    df = groupby(df, [:year, :age, :province, :projection_scenario])
    df = select(df, :, :N => ((x) -> sum(x)) => :N_age, ungroup=false)
    df = transform(df, [:N, :N_age] => ByRow((x, y) -> x / y) => :prop_male, ungroup=true)

    # get the total number of births for a given year, province, and projection scenario
    df_birth = filter(([:age] => (x) -> x == 0), df)
    df_birth = transform(df_birth, [:N_age] => ByRow((x) -> x) => :N_birth)
    df_birth = select(df_birth, Not([:age, :N, :N_age, :prop_male]))

    # add the births column to the main df
    df = leftjoin(df, df_birth, on=[:province, :sex, :year, :projection_scenario])
    df = transform(df, [:N_age, :N_birth] => ByRow((x, y) -> x / y) => :prop)

    # keep only male entries
    df = filter(([:sex] => (x) -> x == "M"), df)
    df = select(df, Not([:sex]))
    return df
end


function process_birth_estimate_data()
    past_population_data = load_past_population_data()
    min_year = findmax(past_population_data.year)[1] + 1
    projected_population_data = load_projected_population_data(min_year)
    birth_estimate = vcat(past_population_data, projected_population_data)
    file_path = joinpath(PROCESSED_DATA_PATH, "master_birth_estimate.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, birth_estimate)
end


function process_initial_population_data()
    past_population_data = load_past_initial_population_data()
    min_year = findmax(past_population_data.year)[1]
    projected_population_data = load_projected_initial_population_data(min_year)
    initial_population = vcat(past_population_data, projected_population_data)
    file_path = joinpath(PROCESSED_DATA_PATH, "master_initial_pop_distribution_prop.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, initial_population)
end


process_birth_estimate_data()
process_initial_population_data()
