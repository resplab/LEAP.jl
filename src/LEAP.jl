module LEAP

# Write your package code here
using DataFrames, Query, CSV, JLD, JLD2, FileIO
using Setfield, Distributions, StatsFuns, StatsBase, Random, SpecialFunctions
using TimerOutputs, Printf
using GRIB

# using Plots
include("global_variables.jl")
include("modules/abstractModule.jl")
include("modules/agent.jl")
include("modules/familyHistory.jl")
include("modules/exacerbation.jl")
include("modules/severity.jl")
include.(filter(contains(r".jl$"), readdir(joinpath(dirname(pathof(LEAP)),"modules/"); join=true)))
include("utils.jl")
include("simulation.jl")


export
    # functions
    process,
    process_initial,
    random_parameter_initialization!,
    create_event_dict,
    set_up,
    create_agent_newborn,
    get_num_newborn,
    agent_has_asthma,
    compute_num_exacerbations,
    # global datasets
    master_birth_estimate,
    master_life_table,
    master_population_initial_distribution,
    master_immigration_table,
    master_emigration_table,
    master_reassessment,
    master_dx,
    master_mis_dx,
    M3_calibrated_asthma_prev_inc
    exacerbation_calibration
    eq5d
end # module
