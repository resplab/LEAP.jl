PROCESSED_DATA_PATH = joinpath(dirname(pathof(LEAP)), "processed_data")

# Antibiotic Exposure
abx_mid_trends = groupby(CSV.read(
    joinpath(PROCESSED_DATA_PATH,"midtrends.csv"),
    DataFrame
),[:year,:sex])

# Asthma Utility
eq5d = CSV.read(joinpath(PROCESSED_DATA_PATH, "eq5d_canada.csv"), DataFrame)
eq5d = groupby(eq5d,[:age,:sex])

# Birth
master_birth_estimate = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_birth_estimate.csv"),
    DataFrame
)
master_population_initial_distribution = CSV.read(
    joinpath(PROCESSED_DATA_PATH,"master_initial_pop_distribution_prop.csv"),
    DataFrame
)

# Census Divisions
master_census_data = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_census_data_2021.csv"),
    DataFrame
)

# Death
master_life_table = CSV.read(joinpath(PROCESSED_DATA_PATH, "master_life_table.csv"), DataFrame)

# Diagnosis
master_dx = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_dx.csv"),
    DataFrame
)
master_mis_dx = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_mis_dx.csv"),
    DataFrame
)

# Emigration
master_emigration_table = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_emigration_table.csv"),
    DataFrame
)

# Exacerbation
exacerbation_calibration  = CSV.read(
    joinpath(PROCESSED_DATA_PATH,"master_calibrated_exac.csv"),
    DataFrame
)

# Immigration
master_immigration_table = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_immigration_table_modified.csv"),
    DataFrame
)

# Incidence
master_incidence_rate = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_inc_interpolated.csv"),
    DataFrame
)
master_prevalence_rate = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_prev_interpolated.csv"),
    DataFrame
)
M3_calibrated_asthma_prev_inc = CSV.read(
    joinpath(PROCESSED_DATA_PATH,"master_calibrated_asthma_prev_inc_M3.csv"),
    DataFrame
)

# Reassessment
master_reassessment = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_assessment.csv"),
    DataFrame
)