using CSV, DataFrames, Logging
include("utils.jl")

PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
STARTING_YEAR = 2000
STARTING_YEAR_PROJ = 2021
MAX_YEAR = 2065
PROVINCES = ["CA", "BC"]



function get_prev_year_population(df, sex, year, age, min_year, min_age)
    if year == min_year || age == min_age
        df = DataFrame(:year => missing, :age => missing, :N => missing, :prob_death => missing)
        return df[1, :]
    else
       return filter(
            ([:sex, :year, :age] => (x, y, z) ->
                    x == sex &&
                    y == year - 1 &&
                    z == age - 1
            ), df
        )[1, [:year, :age, :N, :prob_death]]
    end
end


function get_delta_n(n, n_prev, prob_death)
    return n - n_prev * (1 - prob_death)
end


function get_n_immigrants(n_migrants)
    return n_migrants < 0 ? 0 : n_migrants
end


function load_migration_data()
    @info "Loading initial population data from CSV file..."
    df_population = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_initial_pop_distribution_prop.csv"),
        DataFrame
    )
    @info "Loading mortality data from CSV file..."
    life_table = CSV.read(
        joinpath(PROCESSED_DATA_PATH, "master_life_table.csv"),
        DataFrame
    )

    df_immigration = DataFrame(
        :year => Integer[],
        :province => AbstractString[],
        :age => Integer[],
        :sex => AbstractString[],
        :projection_scenario => AbstractString[],
        :n_immigrants => Float64[],
        :prop_immigrants_birth => Float64[],
        :prop_immigrants_year => Float64[]
    )

    for province in PROVINCES
        @info "Processing data for province $province..."

        # select the required columns
        df = select(filter(
                ([:year, :province] => (x, y) ->
                    x >= STARTING_YEAR &&
                    y == province
                ),
                df_population
            ),
            [:year, :age, :province, :n_age, :prop_male, :projection_scenario]
        )
        df_male = transform(df, [:n_age, :prop_male] => ByRow((x, y) -> x * y) => :N)
        df_male[:, "sex"] = repeat(["M"], size(df_male)[1])
        df_female = transform(df, [:n_age, :prop_male] => ByRow((x, y) -> x * (1 - y)) => :N)
        df_female[:, "sex"] = repeat(["F"], size(df_female)[1])
        df = vcat(df_male, df_female)
        df = select(df, Not([:prop_male, :n_age]))

        projection_scenarios = Set(
            filter(([:projection_scenario] => (x) -> x != "past"), df).projection_scenario
        )
        min_year = findmin(df.year)[1]
        min_age = 0

        for projection_scenario in projection_scenarios
            @info "Projection scenario: $projection_scenario"

            # select only the current projection scenario and the past projection scenario
            df_proj = filter(
                ([:projection_scenario, :year] => (x, y) ->
                x in ["past", projection_scenario] && !(x == "past" && y == STARTING_YEAR_PROJ)
                ),
                df
            )

            # join to the life table to get death probabilities
            df_proj = leftjoin(df_proj, life_table, on=[:year, :age, :province, :sex])

            # get the number of births in each year
            df_birth = filter(([:age] => (x) -> x == 0), df_proj)
            df_birth = groupby(df_birth, [:year])
            df_birth = select(df_birth, :, :N => ((x) -> sum(x)) => :n_birth, ungroup=true)
            df_birth = filter(([:sex] => (x) -> x == "F"), df_birth)
            df_birth = select(df_birth, [:year, :n_birth])

            # get the next year and age for each entry
            df_proj = transform(
                df_proj,
                [:sex, :year, :age] =>
                    ByRow(
                        (x, y, z) -> get_prev_year_population(df_proj, x, y, z, min_year, min_age)
                    )
                    => [:year_prev, :age_prev, :n_prev, :prob_death_prev]
            )

            # remove the missing data
            df_proj = dropmissing(df_proj, :n_prev)

            # compute the change in population
            df_proj = transform(
                df_proj,
                [:N, :n_prev, :prob_death_prev] =>
                    ByRow((x, y, z) -> get_delta_n(x, y, z))
                    => :delta_N
            )

            # add the :n_birth column to df_proj
            df_proj = leftjoin(df_proj, df_birth, on=[:year])

            df_immigration_proj = transform(df_proj, [:delta_N] =>
                ByRow((x) -> get_n_immigrants(x)) => :n_immigrants
            )

            df_immigration_proj = transform(df_immigration_proj, [:n_immigrants, :n_birth] =>
                ByRow((x, y) -> x / y) => :prop_immigrants_birth
            )

            df_immigration_proj = select(
                df_immigration_proj,
                [
                    :province, :year, :age, :sex, :projection_scenario,
                    :prop_immigrants_birth, :n_immigrants
                ]
            )

            # get the immigrants for a given age and sex relative to the immigrants for that year
            df_immigration_proj = groupby(df_immigration_proj, [:year])
            df_immigration_proj = select(
                df_immigration_proj,
                :,
                :n_immigrants => ((x) -> sum(x)) => :n_immigrants_year,
                ungroup=false
            )
            df_immigration_proj = transform(
                df_immigration_proj,
                [:n_immigrants, :n_immigrants_year] =>
                    ByRow((x, y) -> x / y) => :prop_immigrants_year,
                ungroup=true
            )

            df_immigration_proj = select(df_immigration_proj, Not([:n_immigrants_year]))
            df_immigration_proj = transform(
                df_immigration_proj,
                [:projection_scenario] => ByRow((x) -> projection_scenario) => :projection_scenario
            )

            append!(df_immigration, df_immigration_proj)
        end
    end

    file_path = joinpath(PROCESSED_DATA_PATH, "migration/master_immigration_table2.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, df_immigration)
end


load_migration_data()
