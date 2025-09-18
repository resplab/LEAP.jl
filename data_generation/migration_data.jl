using CSV, DataFrames, Logging
include("utils.jl")

PROCESSED_DATA_PATH = joinpath(dirname(dirname(@__FILE__)), "src/processed_data")
DATA_PATH = joinpath(dirname(@__FILE__), "public_dataset")
STARTING_YEAR = 2000
STARTING_YEAR_PROJ = 2021
MAX_YEAR = 2065
PROVINCES = ["CA", "BC"]



function get_prev_year_population(
    df::DataFrame, sex::AbstractString, year::Integer, age::Integer,
    min_year::Integer, min_age::Integer
)
    """
        get_prev_year_population(df, sex, year, age, min_year, min_age)

    Get the age, sex, probability of death, and population for the previous year.

    # Arguments
    - `df::DataFrame`: TODO.
    - `sex::AbstractString`: One of "F" = female, "M" = male.
    - `year::Integer`: The calenar year.
    - `age::Integer`: The integer age.
    - `min_year::Integer`: The minimum year in the dataframe.
    - `min_age::Integer`: The minimum age in the dataframe.

    # Returns
    - `DataFrameRow`: The age, sex, probability of death, and population for the previous year.
    """
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


function get_delta_n(n::Float64, n_prev::Float64, prob_death::Float64)::Float64
    """
        get_delta_n(n, n_prev, prob_death)

    Get the population change due to migration for a given age and sex in a single year.

    # Arguments
    - `n::Float64`: The number of people living in Canada for a single age, sex, and year.
    - `n_prev::Float64`: The number of people living in Canada for the same sex as `n`, in the
        previous year and age. So if `n` is the number of females aged 10 in the year 2020,
        `n_prev` is the number of females aged 9 in the year 2019.
    - `prob_death::Float64`: The probability that a person with a given age and sex in a given
        year will die between the previous year and this year. So if the person is a female
        aged 10 in 2020, `prob_death` is the probability that a female aged 9 in 2019 will die
        by the age of 10.

    # Returns
    - `Float64`: The change in population for a given year, age, and sex due to migration.
    """
    return n - n_prev * (1 - prob_death)
end


function get_n_migrants(delta_N::Float64)::Array{1, Float64}
    """
        get_n_migrants(delta_N)

    Get the number of immigrants and emigrants in a single year for a given age and sex.

    TODO: This function is wrong. delta_N is the change in population due to migration. This
        function currently assumes that if delta_N is less than zero, 100% of migration is
        emigration, and if it is greater than zero, 100% of migration is immigration. This has
        led to the data being very inaccurate (for example, it appears as though people in their
        90s are emigrating a lot and people in their 20s are not). This will be remedied in a
        separate PR.

    # Arguments
    - `delta_N::Float64`: The change in population for a given year, age, and sex due to migration.

    # Returns
    - `Array{1, Float64}`: A vector containing two values, the number of immigrants in a single
        year and the number of emigrants in a single year.
    """
    return [delta_N < 0 ? 0 : delta_N, delta_N > 0 ? 0 : -delta_N]
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

    df_emigration = DataFrame(
        :year => Integer[],
        :province => AbstractString[],
        :age => Integer[],
        :sex => AbstractString[],
        :projection_scenario => AbstractString[],
        :n_emigrants => Float64[],
        :prop_emigrants_birth => Float64[],
        :prop_emigrants_year => Float64[]
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

            # get the number of immigrants/emigrants
            df_migration_proj = transform(df_proj, [:delta_N] =>
                ByRow((x) -> get_n_migrants(x)) => [:n_immigrants, :n_emigrants]
            )

            # compute the proportion of immigrants/emigrants to the number of births in a year
            df_migration_proj = transform(
                df_migration_proj,
                [:n_immigrants, :n_emigrants, :n_birth] =>
                    ByRow(
                        (x, y, z) -> [x / z, y / z]
                    ) => [:prop_immigrants_birth, :prop_emigrants_birth]
            )

            df_migration_proj = select(
                df_migration_proj,
                [
                    :province, :year, :age, :sex, :projection_scenario,
                    :prop_immigrants_birth, :prop_emigrants_birth,
                    :n_immigrants, :n_emigrants
                ]
            )

            # get the migrants for a given age and sex relative to the migrants for that year
            df_migration_proj = groupby(df_migration_proj, [:year])
            df_migration_proj = select(
                df_migration_proj,
                :,
                :n_immigrants => ((x) -> sum(x)) => :n_immigrants_year,
                :n_emigrants => ((x) -> sum(x)) => :n_emigrants_year,
                ungroup=false
            )
            df_migration_proj = transform(
                df_migration_proj,
                [:n_immigrants, :n_immigrants_year, :n_emigrants, :n_emigrants_year] =>
                    ByRow(
                        (x, y, z, w) -> [x / y, z / w]) =>
                        [:prop_immigrants_year, :prop_emigrants_year],
                ungroup=true
            )

            # remove :n_immigrants_year, :n_emigrants_year
            df_migration_proj = select(
                df_migration_proj, Not([:n_immigrants_year, :n_emigrants_year])
            )

            # convert the "past" projection scenario to the given projection scenario
            df_migration_proj = transform(
                df_migration_proj,
                [:projection_scenario] => ByRow((x) -> projection_scenario) => :projection_scenario
            )

            # create separate immigration and emigration dataframes
            df_immigration_proj = select(
                df_migration_proj,
                Not([:n_emigrants, :prop_emigrants_year, :prop_emigrants_birth])
            )
            df_emigration_proj = select(
                df_migration_proj,
                Not([:n_immigrants, :prop_immigrants_year, :prop_immigrants_birth])
            )

            # append the immigration and emigration dataframes for the current projection scenario
            append!(df_immigration, df_immigration_proj)
            append!(df_emigration, df_emigration_proj)
        end
    end

    file_path = joinpath(PROCESSED_DATA_PATH, "migration/master_immigration_table.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, df_immigration)

    file_path = joinpath(PROCESSED_DATA_PATH, "migration/master_emigration_table.csv")
    @info "Saving data to $file_path"
    CSV.write(file_path, df_emigration)
end


load_migration_data()
