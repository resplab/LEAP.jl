PROCESSED_DATA_PATH = joinpath(dirname(pathof(LEAP)), "processed_data")
CONFIG_PATH = joinpath(dirname(pathof(LEAP)), "config.json")

# Antibiotic Exposure
abx_mid_trends = groupby(CSV.read(
    joinpath(PROCESSED_DATA_PATH,"midtrends.csv"),
    DataFrame
),[:year,:sex])


# Diagnosis
master_dx = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_dx.csv"),
    DataFrame
)
master_mis_dx = CSV.read(
    joinpath(PROCESSED_DATA_PATH, "master_asthma_mis_dx.csv"),
    DataFrame
)
