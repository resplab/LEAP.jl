library(tidyverse)
library(mgcv)
library(here)
# goal: estimate the "un"diagnosis prob,
# i.e., the prob that an asthma patient "outgrows" asthma and becomes undiagnosed with asthma

starting_year <- 1999
end_year <- 2065
STABILIZATION_YEAR <- 2025

predict_asthma_occurrence <- function(
    model, age, asthma_max_age, sex, year, stabilization_year=STABILIZATION_YEAR
) {
    age <- pmin(age, asthma_max_age)
    year <- pmin(year, stabilization_year)
    return(exp(predict(model, newdata=data.frame(age, sex, year))) %>% unlist())
}


get_reassessment_data <- function(
    chosen_province="CA",
    starting_year=1999,
    end_year=2065,
    stabilization_year=2025,
    asthma_max_age=62
) {

    asthma_inc_model <- read_rds(here("R/asthma_incidence_model.rds"))
    asthma_prev_model <- read_rds(here("R/asthma_prevalence_model.rds"))

    df_asthma <- expand.grid(age=3:110, sex=c(0, 1), year=starting_year:end_year) %>%
        as.data.frame()

    df_asthma <- df_asthma %>%
        mutate(
            inc=predict_asthma_occurrence(
                asthma_inc_model, age, asthma_max_age, sex, year,
                stabilization_year
            ),
            prev=predict_asthma_occurrence(
                asthma_prev_model, age, asthma_max_age, sex, year,
                stabilization_year
            )
        ) %>%
        # let inc be prev for age = 3
        mutate(inc=ifelse(age == 3, prev, inc))

    df_asthma_year <- df_asthma %>% group_split(year)

    results_assessment <- c()
    years <- c((starting_year + 1):end_year)
    ages <- c(4:110)

    for (i in 1:(length(years) - 1)) {

        tmp_year <- years[i]

        # Get the predicted prevalence for the previous year
        df_prevalence_past <- df_asthma_year[[i]] %>%
            select(age, year, sex, prev) %>%
            filter(age != 110) %>%
            pivot_wider(names_from=sex, values_from=prev) %>%
            select(-age, -year)

        # Get the predicted prevalence for the current year
        df_prevalence <- df_asthma_year[[i + 1]] %>%
            select(age, year, sex, prev) %>%
            filter(age != 3) %>%
            pivot_wider(names_from=sex, values_from=prev) %>%
            select(-age, -year)

        # Get the predicted incidence for the current year
        df_incidence <- df_asthma_year[[i + 1]] %>%
            select(age, year, sex, inc) %>%
            filter(age != 3) %>%
            pivot_wider(names_from=sex, values_from=inc) %>%
            select(-age, -year)

        # estimate tmp_u by assuming P(correct diagnosis) = 1
        df_assessment <- (
            df_prevalence - df_incidence * (1 - df_prevalence_past)
        ) / df_prevalence_past

        results_assessment[[i]] <- cbind(
            data.frame(age=4:110, year=year),
            df_assessment
        )
    }

    df_reassessment <- results_assessment %>%
        do.call(rbind, .) %>%
        as.data.frame()

    colnames(df_reassessment)[c(3, 4)] <- c("F", "M")
    df_reassessment$province <- chosen_province

    return(df_reassessment)
}


provinces <- c("BC", "CA")
end_years <- c(2043, 2066)
reassessment_list <- c()
for (i in length(provinces)) {
    df <- get_reassessment_data(
        province=provinces[i], end_year=end_years[i], stabilization_year=STABILIZATION_YEAR
    )
    reassessment_list[[i]] <- df
}

df_reassessment <- do.call(rbind, reassessment_list) %>% as.data.frame()

df_reassessment <- df_reassessment %>% pivot_longer(3:4, names_to="sex", values_to="reassessment")
df_reassessment <- df_reassessment %>%
    mutate(reassessment=ifelse(reassessment > 1, 1, reassessment))

write_csv(df_reassessment, here("src/processed_data/asthma_reassessment.csv"))
