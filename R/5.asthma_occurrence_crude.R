# source: https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1310009608

library(tidyverse)
library(here)
library(mgcv)
max_year <- 2019
STARTING_YEAR <- 2000
STABILIZATION_YEAR <- 2025
MIN_AGE <- 3


#' Load the asthma incidence and prevalence data from BC HSDA dataset
#'
#' Data Columns:
#' agegrp: str, format "X-Y", "80+"
#' Year: int, format XXXX
#' SEX: str, "M", "F", or "T"
#' incidence
#' prevalence
#' incidence_numerator
#' prevalence_numerator
#' pop: int, number of people in the age group
#' @param starting_year An integer. The starting year of the data. Default is 2000.
#' @param min_age An integer. The minimum age for asthma diagnosis. Default is 3.
#' @return Returns an object of class "?". Description of what the function returns
#' @details This function can be used in place of the `load_asthma_df_bc` function if needed.
load_asthma_df_hsda <- function(starting_year=STARTING_YEAR, min_age=MIN_AGE) {
    df <- read_csv(here("R/private_dataset/asthma_incidence_prevalence_bc_hsda_2001_2023.csv"))
    # Rename columns
    df <- df %>% rename(age_group=agegrp)

    # Filter for year >= starting_year
    df <- df %>% filter(year >= starting_year)

    # Age groups are in the format "X-Y" or "X+"
    # Set the age to the average of the age group
    lapply(df$age_group, function(x){
        ceiling(mean(parse_number(str_split(x, "-")[[1]])))
    })  %>% 
    unlist() -> df$age

    # Key assumption: asthma starts at age 3
    # Set incidence = prevalence at age 3
    df <- df %>% 
        mutate(
            incidence=ifelse(age==min_age, PREV_CRUDE_RATE, INC_CRUDE_RATE)
        )

    # AH
    df <- df %>% 
        filter(REGION == "BC") %>%
        filter(DISEASE == "Asthma 1+") %>%
        filter(SEX != "Total") %>%
        filter(age_group != "Total") %>%
        rename(year=Year, sex=SEX)

    return(df)
}


#' Load the asthma incidence and prevalence data from BC administrative dataset
#'
#' Data Columns:
#' age_group_desc: str, format "X-Y years", "<1 year", "90+ years"
#' fiscal_year: str, format XXXX[A-z0-9]
#' gender: str, "M", "F", or "T"
#' incidence
#' prevalence
#' incidence_numerator
#' prevalence_numerator
#' pop: int, number of people in the age group
#' @param starting_year An integer. The starting year of the data. Default is 2000.
#' @return Returns an object of class "?". Description of what the function returns
load_asthma_df_bc <- function(starting_year=STARTING_YEAR) {
    df <- readxl::read_xlsx(here("R/private_dataset/asthma_inc_prev.xlsx"), sheet=1)

    # Rename columns
    df <- df %>% rename(age_group=age_group_desc, sex=gender)

    # Filter out the "<1 year" age group
    df <- df %>%
        filter(age_group != "<1 year")

    df <- df %>% 
        mutate(year=substr(fiscal_year, 1, 4) %>% as.numeric())

    # Filter for year >= starting_year
    df <- df %>%
        filter(year >= starting_year)

    df <- df %>%
        mutate(
            prev_sd=sqrt(prevalence_numerator) * qnorm(0.975),
            prev_upper=(prevalence_numerator + prev_sd)/pop,
            prev_lower=(prevalence_numerator - prev_sd)/pop
        )

    # Age groups are in the format "X-Y years" or "<1 year"
    # Set the age to the average of the age group
    lapply(df$age_group,function(x){
        ceiling(mean(parse_number(str_split(x, "-")[[1]])))
    })  %>% 
    unlist() -> df$age

    # Key assumption: asthma starts at age 3
    # Set incidence = prevalence at age 3
    df <- df %>% 
        mutate(
            incidence=ifelse(age==3, prevalence, incidence),
            incidence_numerator=ifelse(age==3, prevalence_numerator, incidence_numerator)
        )

    return(df)
}


plot_occurrence_comparison <- function(
    df, df_pred, title, year_min=2000, year_max=2020, year_step=2
) {
    years <- seq(year_min, year_max, by=year_step)
    ggplot(
        data=df %>% 
            filter(year %in% years) %>% 
            drop_na() %>%
            mutate(year=as.factor(year)),
        aes(x=age, y=y, color=year)
    ) +
        geom_line() +
        geom_line(
            data=df_pred %>% 
                filter(year %in% years) %>% 
                mutate(year=as.factor(year)),
            aes(x=age, y=y, color=year),
            linetype="dashed"
        ) +
        ylab(title) +
        facet_grid(.~sex)
}


#' Generate a generalized additive model (GAM) for asthma incidence in BC
#' 
#' Formula: log(incidence) ~ sex(year + sex*poly(age, 5))
#'
#' @param df_asthma A data frame containing asthma data with columns:
#' year: int, the calendar year, a value between 2000 and 2020
#' sex: str, "M" or "F"
#' age: int, the age of the individual, a value between 3 and 65
#' incidence: numeric, the incidence rate of asthma per 100 individuals
#' @param min_age An integer. The minimum age for asthma diagnosis. Default is 3.
#' @param max_age An integer. The maximum age for the regression analysis. Default is 65.
#' @param starting_year An integer. The starting year of the data. Default is 2000.
#' @details This function generates a GAM for asthma incidence in BC, plots the model predictions,
#' and saves the model to an RDS file.
generate_incidence_model <- function(
    df_asthma,
    min_age=MIN_AGE,
    max_age=65,
    starting_year=STARTING_YEAR
) {
    # Create incidence dataframe
    df <- df_asthma %>% select(year, sex, age, incidence)

    # Filter for year <= max_year
    df <- df %>%
        filter(year >= starting_year & year <= max_year)

    df <- df %>%
        mutate(sex=as.numeric(sex=="M"), y=incidence)
        filter(age <= max_age)

    model <- mgcv::gam(
        formula=log(incidence)~ sex*year + sex*poly(age, degree = 5), data=df
    )
    print(summary(model))

    df_pred <- expand.grid(year=c(2000:2065), sex=c(0,1), age=seq(min_age, max_age, by=1)) %>% 
        as.data.frame()

    df_pred$y <- exp(predict(model, newdata=df_pred))

    df_inc_pred <- df %>% 
        left_join(df_pred, by=c("year","sex","age"))

    plot_occurrence_comparison(
        df=df,
        df_pred=df_inc_pred, 
        title="Asthma Incidence per 100 in BC",
        year_min=2000,
        year_max=2020
    )

    write_rds(model, here("R/private_dataset/asthma_incidence_model.rds"))
}


#' Generate a generalized additive model (GAM) for asthma prevalence in BC
#'
#' Formula: log(prevalalence) ~ sex * poly(year, degree=2) * poly(age, degree=5)
#' 
#' @param df_asthma A data frame containing asthma data with columns:
#' year: int, the calendar year, a value between 2000 and 2020
#' sex: str, "M" or "F"
#' age: int, the age of the individual, a value between 3 and 65
#' prevalence: numeric, the prevalence rate of asthma per 100 individuals
#' @param min_age An integer. The minimum age for asthma diagnosis. Default is 3.
#' @param max_age An integer. The maximum age for the regression analysis. Default is 65.
#' @param starting_year An integer. The starting year of the data. Default is 2000.
#' @details This function generates a GAM for asthma prevalence in BC, plots the model predictions,
#' and saves the model to an RDS file.
generate_prevalence_model <- function(
    df_asthma,
    min_age=MIN_AGE,
    max_age=65,
    starting_year=STARTING_YEAR
) {
    # Create prevalence dataframe
    df <- df_asthma %>% select(year, sex, age, prevalence)

    # Filter for year <= max_year
    df <- df %>% filter(year <= max_year) %>% 
        filter(year >= starting_year & year <= max_year)
    df <- df %>%
        mutate(
            sex=as.numeric(sex=="M"),
            y=prevalence
        ) %>%
        filter(age <= max_age)

    model <- mgcv::gam(
        formula=log(prevalence)~ sex*poly(year, degree=2)*poly(age, degree=5), data=df
    )

    print(summary(model))

    df_pred <- expand.grid(year=c(2000:2065), sex=c(0, 1), age=seq(min_age, 62, by=1)) 

    df_pred$y <- exp(predict(model, newdata=df_pred))

    df_prev_pred <- df %>% 
        left_join(df_pred, by=c("year","sex","age"))

    plot_occurrence_comparison(
        df=df,
        df_pred=df_prev_pred, 
        title="Asthma Prevalence per 100 in BC",
        year_min=2000,
        year_max=2025
    )

    write_rds(model, here("R/private_dataset/asthma_prevalence_model.rds"))
}


df_asthma <- load_asthma_df_bc()
generate_incidence_model(df_asthma=df_asthma, min_age=MIN_AGE, max_age=65)
generate_prevalence_model(df_asthma=df_asthma, min_age=MIN_AGE, max_age=65)




# generate crude asthma inc & prev  ----------------------------------------

prev_model <- read_rds(here("R/asthma_prevalence_model.rds"))
inc_model <- read_rds(here("R/asthma_incidence_model.rds"))
max_age <- 63
df <- expand.grid(year=STARTING_YEAR:2065, sex=c(0:1), age=MIN_AGE:110) %>% 
    as.data.frame() %>% 
    mutate(
        prev=exp(predict(
            prev_model,
            data.frame(year=pmin(STABILIZATION_YEAR, year), sex, age=pmin(age, max_age))
        )),
        inc=exp(predict(
            inc_model,
            data.frame(year=pmin(STABILIZATION_YEAR, year), sex, age=pmin(age, max_age))
        ))
    ) %>% 
    mutate(
        prev=as.numeric(prev),
        inc=as.numeric(inc)
    )

write_csv(df, here("R/master_asthma_prev_inc.csv"))
ggplot(
    data=df %>% 
        mutate(sex = ifelse(sex==1, "Male", "Female")) %>% 
        filter(year %in% seq(STARTING_YEAR, STABILIZATION_YEAR, by=5)) %>% 
        mutate(year=as.factor(year)),
    aes(x=age, y=prev, col=year)
) +
    geom_line() +
    xlim(c(0, 60)) +
    facet_grid(.~sex) +
    ylab("Crude asthma prevalence (per 100)")+
    xlab("Age (years)")+
    theme_classic(base_size=20) +
    theme(
        legend.position='top',
        legend.title=element_blank()
    )

ggplot(
    data=df %>% 
        mutate(sex=ifelse(sex==1, "Male", "Female")) %>% 
        filter(year %in% seq(STARTING_YEAR, STABILIZATION_YEAR, by=5)) %>% 
        mutate(year=as.factor(year)),
    aes(x=age, y=inc, col=year)
) +
    geom_line() +
    xlim(c(0, 60)) +
    facet_grid(.~sex) +
    ylab("Crude asthma incidence (per 100)")+
    xlab("Age (year)")+
    theme_classic(base_size=20) +
    theme(
        legend.position='top',
        legend.title=element_blank()
    )

write_csv(df, here("R/master_asthma_prev_inc.csv"))

