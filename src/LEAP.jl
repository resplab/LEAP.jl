module LEAP

# Write your package code here
using DataFrames, Query, CSV, JLD, JLD2, FileIO, JSON
using Setfield, Distributions, StatsFuns, StatsBase, Random, SpecialFunctions
using TimerOutputs, Printf
using GRIB, PolygonOps, Shapefile, JSON, GeoInterface

PROCESSED_DATA_PATH = joinpath(dirname(pathof(LEAP)), "processed_data")
CONFIG_PATH = joinpath(dirname(pathof(LEAP)), "config.json")

include("modules/abstractmodule.jl")
include("modules/birth.jl")
include("modules/familyhistory.jl")
include("modules/antibioticexposure.jl")
include("modules/agent.jl")
include("modules/exacerbation.jl")
include("modules/severity.jl")
include.(filter(contains(r".jl$"), readdir(joinpath(dirname(pathof(LEAP)),"modules/"); join=true)))
include("utils.jl")
include("outcomematrix.jl")
include("simulation.jl")


export
    # functions
    create_agent,
    process,
    process_initial,
    random_parameter_initialization!,
    create_outcome_matrix,
    set_up,
    get_num_newborn,
    get_initial_population_indices,
    get_num_new_immigrants,
    get_new_agents,
    agent_has_asthma,
    compute_num_exacerbations,
    load_census_boundaries,
    point_in_polygon,
    # global datasets
    PROCESSED_DATA_PATH,
    CONFIG_PATH
end # module
