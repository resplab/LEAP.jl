using GRIB
using DataFrames
using Statistics
using CSV


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
    year::Integer
    month::Union{Nothing, Integer}
    day::Union{Nothing, Integer}
    projection::String
    longitudes::Array{Float64, 1}
    latitudes::Array{Float64, 1}
    values::Array{Float64, 1}
end


"""
    load_grib_file(filename)

Load a *.grib2 file and almalgamate the data by taking the mean.

# Arguments
- `filename::AbstractChar`: the .grib2 file to open.

"""
function load_grib_file(filename::AbstractChar)
    df = DataFrame()
    GribFile(filename) do f
        index = 1
        for record in f
            df = add_record_to_df!(df, record["longitudes"], record["latitudes"],
                record["codedValues"], index)
            index += 1
            year = record["year"]
            month = record["month"]
            day = record["day"]
            projection = record["gridType"]
        end
    end

    grib_data = get_grib_data_average(df, year, month, day, projection)
    return grib_data
end


"""
    load_grib_files(folder)

Load multiple *.grib2 files and almalgamate the data by taking the mean.

# Arguments
- `folder::AbstractChar`: the folder containing the .grib2 files to open.

"""
function load_grib_files(folder::AbstractChar)
    index = 1
    df = DataFrame()
    for filename in readdir(folder)
        extension = filename[findlast(isequal('.'), filename):end]
        if extension == ".grib2"
            grib_data = load_grib_file(filename)
            df = add_record_to_df!(df, grib_data.longitudes, grib_data.latitudes,
                grib_data.values, index)
            index += 1
        end
    end
    final_grib_data = get_grib_data_average(df, year, month, day, projection)
    return final_grib_data
end



"""
    add_record_to_df(df, longitudes, latitudes, values, index)

Add a new column to the grib data frame.

# Arguments
- `df::DataFrame`: a data frame with the following columns:
    `longitudes`: Float64, a list of longitude values.
    `latitudes`: Float64, a list of latitude values.
    `value_$index`: Float64, a list of the values of interest at a specified longitude and latitude.
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
    `value_$index`: Float64, a list of the values of interest at a specified longitude and latitude.
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


function save_grib_data(grib_data::GribData, file_path::String)
    df = DataFrame()
    df[!, "longitudes"] = grib_data.longitudes
    df[!, "latitudes"] = grib_data.latitudes
    df[!, "values"] = grib_data.values
    CSV.write(file_path, df)
end
