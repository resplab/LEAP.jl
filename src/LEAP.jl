module LEAP

# Write your package code here
using DataFrames, Query, CSV, JLD2, FileIO
using Setfield, Distributions, StatsFuns, StatsBase, Random, SpecialFunctions
using TimerOutputs, Printf

# using Plots

include("simulation.jl")

export
    # functions
    process,
    process_initial,
    random_parameter_initialization!,
    set_up,
    # global datasets
    birth_projection,
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
