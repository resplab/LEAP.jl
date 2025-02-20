library(tidyverse)
library(here)
library(mgcv)
library(roptim)

# we need: 
# 1) asthma prev
# 2) target rate (this is in per pop)
# 3) population 
# 4) 2) + 3) => annual number of severe hospitalizations needed
# 5) 1) + exacerbation module => annual number of exacerbations
# 6) so we can identify a tuner to do this

# asthma prev


logistic <- function(x){
    if(is.infinite(exp(x))){
        return(1)
    } else {
        return(exp(x)/(1+exp(x)))
    }
}

control_prediction <- function(sex, age, theta=c(-1e5,-0.3950, 2.754,1e5)){
    age_scaled = age / 100
    eta <- age_scaled*3.5430381 + sex*0.2347807 + age_scaled * sex *-0.8161495 +
        age_scaled^2 * sex *-1.1654264 + age_scaled^2 * -3.4980710
    results <- rep(0,3)
    for(j in 1:(length(theta)-1)){
        results[j] <- logistic(theta[j+1]-eta) - logistic(theta[j]-eta)
    }
    results
}

exacerbation_prediction <- function(
    sex, age, beta_control=log(c(0.1880058, 0.3760116, 0.5640174))
){
    if(age<3){
        return(0)
    }
    control <- control_prediction(sex, age)
    exp(sum(control * beta_control))
}


exacerbation_calibrator <- function(
    chosen_province="CA", baseline_year=2000, max_cal_year=2065, stablization_year=2025
){

    min_cal_year <- baseline_year

    if(chosen_province=="CA"){
        max_cal_year <- 2065
    } else {
        max_cal_year <- 2043
    }
    chosen_projection_scenario <- "M3"
    growth_type <- "M3"

    master_prev_inc <- read_csv(here("src/processed_data/master_asthma_prev_inc.csv"))
    df_prev <- master_prev_inc %>%
        select(year, age, sex, prev)

    df_inc <- master_prev_inc %>%
        select(year, sex, age, inc)

    # target population
    df_cihi <- read_rds(
        here(paste0("src/processed_data/asthma_hosp/", chosen_province, "/tab1.rds"))
    )$rate
    df_cihi <- df_cihi %>% 
        filter(fiscal_year >= baseline_year) %>% 
        rename(year=fiscal_year)
    df_cihi <- df_cihi %>% 
        pivot_longer(-1, names_to="type", values_to="true_rate")
        
    df_cihi <- df_cihi %>% 
        filter(!is.na(true_rate)) 
        
    # Remove the "+" from the type entries
    df_cihi <- df_cihi %>% 
        mutate(type=str_remove(type,"\\+"))
        
    df_cihi <- df_cihi %>% 
        mutate(sex=case_when(
                str_detect(type,"M") ~ 1,
                str_detect(type,"F") ~ 0,
                TRUE ~ NA)
        )
    df_cihi <- df_cihi %>% 
        filter(!is.na(sex))
        
    df_cihi <- df_cihi %>% 
        mutate(age=parse_number(type))
        
    df_cihi <- df_cihi %>% 
        filter(!is.na(age)) %>% 
        select(-type) %>% 
        filter(age>=3) %>% 
        arrange(year,sex,age)

    # Get the most recent year df
    tmp_cihi <- df_cihi %>% filter(year == max(df_cihi$year))

    impute_years <- (max(df_cihi$year)+1):(max_cal_year)

    for(i in impute_years){
        df_cihi <- rbind(df_cihi, tmp_cihi %>% mutate(year=i))
    }

    # past population
    df_population <- read_csv(here("src/processed_data/master_initial_pop_distribution_prop.csv")) %>% 
        filter(province==chosen_province) %>% 
        filter(year >= baseline_year) %>%
        mutate(sex="M") %>%
        mutate(n_age_sex=n_age * prop_male)
    df_population <- rbind(df_population, df_population %>% 
        mutate(sex="F", n_age_sex=n_age * (1-prop_male)))
    df_population <- df_population %>%
        filter(
            (projection_scenario == "past" & year <= 2021) |
            (projection_scenario == growth_type & year > 2021)
        ) %>%
        select(-projection_scenario, -n_age, -prop_male, -prop, -n_birth) %>%
        rename(n=n_age_sex) %>%
        mutate(case_when(sex=="M" ~ 1, sex=="F" ~ 0))

    pop <- df_population %>% 
        filter(year <= max_cal_year) %>% 
        filter(age>=3)
    pop <- pop %>%
        mutate(age=ifelse(age>90, 90, age))
    pop <- pop %>%
        group_by(year, sex, age) %>% 
        summarise(n=sum(n)) %>%
        ungroup()
    pop <- pop %>%
        mutate(sex=as.numeric(sex=="M"))

    p_hosp <- 0.026

    df_target <- pop %>% 
        left_join(df_cihi, by=c("year", "sex", "age"))

    df_target <- df_target %>%
        mutate(true_n=true_rate * n / 100000)

    df_target <- df_target %>%
        left_join(df_prev, by=c("year","sex","age"))

    df_target <- df_target %>%
        mutate(n_asthma=prev*n)

    df_target <- df_target %>%
        rowwise() %>%
        mutate(
            mean_annual_exacerbation=exacerbation_prediction(sex, age),
            expected_exacerbations=mean_annual_exacerbation*n_asthma,
            expected_n=p_hosp*expected_exacerbations,
            calibrator_multiplier=true_n/expected_n
        )

    exac_cal <- df_target %>% 
        select(year, sex, age, calibrator_multiplier) %>% 
        mutate(province=chosen_province)

    exac_cal
}

exacerbation_calibration_BC <- exacerbation_calibrator("BC")
exacerbation_calibration_CA <- exacerbation_calibrator("CA")

final_result <- rbind(exacerbation_calibration_BC,exacerbation_calibration_CA)

write_csv(final_result, here("src/processsed_data/master_calibrated_exac.csv"))


# Convert the rds files to csv

PROVINCES = c("CA", "AB", "BC", "MB", "NB", "NL", "NS", "ON", "PE", "QC", "SK", "TR")
for(province in PROVINCES){
    print(province)
    df_cihi <- read_rds(
        here(paste0("src/processed_data/asthma_hosp/", province, "/tab1.rds"))
    )
    for(name in names(df_cihi)){
        print(name)
        write_csv(
            df_cihi[[name]],
            here(paste0("src/processed_data/asthma_hosp/", province, "/tab1_", name, ".csv"))
        )
    }
}