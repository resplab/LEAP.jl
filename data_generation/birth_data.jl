using CSV, DataFrames, Logging

DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
STARTING_YEAR = 1999

PROVINCE_MAP = Dict(
    "Canada" => "CA",
    "British Columbia" => "BC",
    "Alberta" => "AB",
    "Saskatchewan" => "SK",
    "Manitoba" => "MB",
    "Ontario" => "ON",
    "Quebec" => "QC",
    "Newfoundland and Labrador" => "NL",
    "Nova Scotia" => "NS",
    "New Brunswick" => "NB",
    "Prince Edward Island" => "PE",
    "Yukon" => "YT",
    "Northwest Territories" => "NT",
    "Nunavut" => "NU"
)

function get_province_id(province)
    return PROVINCE_MAP[province]
end

function get_sex_id(sex)
    return sex[1:1]
end

function get_projection_scenario_id(projection_scenario)
    projection_scenario_id = replace(projection_scenario, "Projection scenario " => "")
    return projection_scenario_id[1:2]
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


function process_birth_estimate_data()
    past_population_data = load_past_population_data()
    min_year = findmax(past_population_data.year)[1] + 1
    projected_population_data = load_projected_population_data(min_year)
    birth_estimate = vcat(past_population_data, projected_population_data)
    file_path = joinpath(PROCESSED_DATA_PATH, "master_birth_estimate.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, birth_estimate)
end

process_birth_estimate_data()
