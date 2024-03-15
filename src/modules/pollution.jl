using GRIB
using DataFrames
using Statistics
using CSV


"""
    PollutionTable

A struct containing information about PM2.5 pollution in Canada.

# Fields
- `data::Union{GroupedDataFrame, Nothing}`: A data frame grouped by the SSP scenario, with the
    following columns:
    `CDUID`: the census division identifier.
    `year`: the year for the pollution data projection.
    `month`: the month for the pollution data projection.
    `date`: the data for the pollution data projection, first of the month.
    `background_pm25`: the average background PM2.5 levels for a given month.
    `wildfire_pm25`: the average PM2.5 levels due to wildfires for a given month.
    `factor`: the future climate scaling factor, based on the SSP scenario.
    `wildfire_pm25_scaled`: `wildfire_pm25` * `factor`.
    `total_pm25`: the total average PM2.5 levels for a given month:
      `wildfire_pm25_scaled` + `background_pm25`
    `SSP`: the SSP scenario, one of `SSP1_2.6`, `SSP2_4.5`, `SSP3_7.0`, `SSP5_8.5`.
"""
@kwdef struct PollutionTable <: PollutionTableModule
    data::Union{GroupedDataFrame{DataFrame}, Nothing}
end


"""
    Pollution

A struct containing information about PM2.5 pollution for a given census division, date, and
SSP scenario.

# Fields
- `cduid::Union{Integer, Nothing}`: the census division identifier.
- `year`: the year for the pollution data projection.
- `month`: the month for the pollution data projection.
- `wildfire_pm25_scaled`: `wildfire_pm25` * `factor`.
- `total_pm25`: the total average PM2.5 levels for a given month:
  `wildfire_pm25_scaled` + `background_pm25`
- `SSP`: the SSP scenario, one of `SSP1_2.6`, `SSP2_4.5`, `SSP3_7.0`, `SSP5_8.5`.
"""
@kwdef struct Pollution <: PollutionModule
    cduid::Union{Integer, Nothing}
    year::Union{Integer, Nothing}
    wildfire_pm25_scaled_mean::Union{Float64, Nothing}
    total_pm25_mean::Union{Float64, Nothing}
    SSP::Union{String, Nothing}
end


"""
    load_pollution_table(pm25_data_path)

Load the data from the PM2.5 SSP *.csv files.

# Arguments
- `pm25_data_path::String`: full directory path for the PM2.5 *.csv files.

# Returns
- `PollutionTable`: an object containing the PM2.5 pollution data for various SSP scenarios.
"""
function load_pollution_table(pm25_data_path::String)
    files = readdir(pm25_data_path)
    pollution_data = DataFrame()
    for file in files
        if splitext(file)[2] == ".csv"
            df = CSV.read(joinpath(pm25_data_path, file), DataFrame)
            pollution_data = [pollution_data;df]
        end
    end
    pollution_table = PollutionTable(
        data = groupby(pollution_data, :SSP)
    )
    return pollution_table
end


"""
    assign_pollution(cduid, year, SSP, pollution_table)

Get the pollution data for a specific year and SSP scenario.

# Arguments
- `cduid::Integer`: the census division identifier.
- `year::Integer`: the year for the pollution data projection.
- `SSP::String`: the SSP scenario, one of `SSP1_2.6`, `SSP2_4.5`, `SSP3_7.0`, `SSP5_8.5`.
- `pollution_table::PollutionTable`: an object containing the PM2.5 pollution data for various
    SSP scenarios.

# Returns
- `Pollution`: an object containing the PM2.5 pollution data for a specific year and SSP scenario.
"""
function assign_pollution(cduid::Integer, year::Integer, SSP::String,
    pollution_table::PollutionTable)

    df = filter(
        [:CDUID, :year] => (x, y) -> x == cduid && y == year, pollution_table.data[(SSP,)]
    )

    pollution = Pollution(
        cduid=cduid,
        year=year,
        wildfire_pm25_scaled_mean=mean(df.wildfire_pm25_scaled),
        total_pm25_mean=mean(df.total_pm25),
        SSP=SSP
    )
    return pollution
end


"""
    GribData

A struct containing GRIB data on air pollution.

# Fields
- `year::Integer`: year the data was collected.
- `month::Integer`: month the data was collected.
- `day::Integer`: day the data was collected.
- `projection::String`: which map projection was used.
    See `gridType` on the [GRIB keys](https://confluence.ecmwf.int/display/ECC/GRIB+Keys) page.
- `longitudes::Array{Float64, 1}`: a list of longitude values.
- `latitudes::Array{Float64, 1}`: a list of latitude values.
- `values::Array{Float64, 1}`: a list of the values of interest at a specified longitude and
    latitude. For example, it could be the PM2.5 concentration.
"""
struct GribData
    year::Union{Nothing, Integer}
    month::Union{Nothing, Integer}
    day::Union{Nothing, Integer}
    projection::Union{Nothing, String}
    longitudes::Union{Nothing, Array{Float64, 1}}
    latitudes::Union{Nothing, Array{Float64, 1}}
    values::Union{Nothing, Array{Float64, 1}}
end


"""
    load_grib_file(filename)

Load a *.grib2 file and almalgamate the data by taking the mean.

# Arguments
- `filename::String`: the .grib2 file to open.

"""
function load_grib_file(filename::String)
    df = DataFrame()
    year = 0
    month = 0
    day = 0
    projection = ""
    GribFile(filename) do f
        index = 1
        for record in f
            longitudes = invert_longitude.(record["longitudes"])
            df = add_record_to_df!(df, longitudes, record["latitudes"],
                record["codedValues"], index)
            index += 1
            year = record["year"]
            month = record["month"]
            day = record["day"]
            projection = String(record["gridType"])
        end
    end
    grib_data = get_grib_data_average(df, year, month, day, projection)
    return grib_data
end


"""
    load_grib_files(folder)

Load multiple *.grib2 files and almalgamate the data by taking the mean.

# Arguments
- `folder::String`: the folder containing the .grib2 files to open.
- `recursive::Bool`: if true, iterate through all subdirectories and compute an
    aggregate average. If false, only read `.grib2` files in the given directory.

"""
function load_grib_files(folder::String, recursive::Bool=false)
    index = 1
    df = DataFrame()
    grib_data = GribData(nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    for item in readdir(folder, join=true)
        if isdir(item)
            if recursive
                println("Reading directory ", item)
                grib_data = load_grib_files(item)
                df = add_record_to_df!(df, grib_data.longitudes, grib_data.latitudes,
                    grib_data.values, index)
                index += 1
            end
        else
            println("Reading file ", item)
            filename = item
            extension = filename[findlast(isequal('.'), filename):end]
            if extension == ".grib2"
                grib_data = load_grib_file(string(filename))
                df = add_record_to_df!(df, grib_data.longitudes, grib_data.latitudes,
                    grib_data.values, index)
                index += 1
            end
        end
    end
    final_grib_data = get_grib_data_average(df, grib_data.year, grib_data.month, nothing,
        grib_data.projection)
    return final_grib_data
end



"""
    add_record_to_df(df, longitudes, latitudes, values, index)

Add a new column to the grib data frame.

# Arguments
- `df::DataFrame`: a data frame with the following columns:
    `longitudes`: Float64, a list of longitude values.
    `latitudes`: Float64, a list of latitude values.
    `value_{index}`: Float64, a list of the values of interest at a specified longitude and latitude.
        For example, it could be the PM2.5 concentration. Each value column corresponds to either
        a record in the original grib file, or a file in a folder.
- `longitudes::Array{Float64, 1}`: a list of longitude values.
- `latitudes::Array{Float64, 1}`: a list of latitude values.
- `values::Array{Float64, 1}`: a list of the values of interest at a specified longitude and
    latitude. For example, it could be the PM2.5 concentration.
- `index::Integer`: the index of the iteration loop.

"""
function add_record_to_df!(df::DataFrame, longitudes::Array{Float64, 1},
    latitudes::Array{Float64, 1}, values::Array{Float64, 1}, index::Integer)
    colname = "value_$index"
    colname_lat = "latitudes"
    colname_lon = "longitudes"
    if index == 1
        df[!, colname_lon] = longitudes
        df[!, colname_lat] = latitudes
        df[!, colname] = values
    else
        df[!, colname] = values
    end
    return df
end


"""
    get_grib_data_average(df, year, month, day, projection)

Find the mean of all the `value_*` columns, and return a `GribData` object.

# Arguments
- `df::DataFrame`: a data frame with the following columns:
    `longitudes`: Float64, a list of longitude values.
    `latitudes`: Float64, a list of latitude values.
    `value_{index}`: Float64, a list of the values of interest at a specified longitude and latitude.
        For example, it could be the PM2.5 concentration. Each value column corresponds to either
        a record in the original grib file, or a file in a folder.
- `year::Integer`: year the data was collected.
- `month::Integer`: month the data was collected.
- `day::Integer`: day the data was collected.
- `projection::String`: which map projection was used.
    See `gridType` on the [GRIB keys](https://confluence.ecmwf.int/display/ECC/GRIB+Keys) page.
"""
function get_grib_data_average(df::DataFrame, year::Integer, month::Union{Nothing, Integer},
    day::Union{Nothing, Integer}, projection::String)
    value_columns = names(select(df, r"value"))
    df = transform(df, AsTable(value_columns) => ByRow(mean) => :mean)
    grib_data = GribData(
        year,
        month,
        day,
        projection,
        df.longitudes,
        df.latitudes,
        df.mean
    )
    return grib_data
end


"""
    save_grib_data(grib_data, file_path)

Given a GribData object, save it to a *.csv file.

# Arguments
- `grib_data::GribData`: a GribData object, see [`GribData`](@ref).
- `file_path::String`: full file name of *.csv file to save the data to.

"""
function save_grib_data(grib_data::GribData, file_path::String)
    df = DataFrame()
    df[!, "longitudes"] = grib_data.longitudes
    df[!, "latitudes"] = grib_data.latitudes
    df[!, "values"] = grib_data.values
    CSV.write(file_path, df)
end


"""
    invert_longitude(longitude)

Longitude values can be given as either + degrees East of the Prime Meridian, or as
- degrees West of the Prime Meridian; convert between the two.

# Arguments
- `longitude::Float64`: a number in [0, 360) giving the degrees east or west of the Prime Meridian.

"""
function invert_longitude(longitude::Float64)
    if longitude > 0
        return longitude - 360
    elseif longitude < 0
        return 360 + longitude
    else
        return longitude
    end
end
