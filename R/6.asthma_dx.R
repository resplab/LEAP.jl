library(tidyverse)
library(mgcv)
# goal: estimate the "un"diagnosis prob, 
# i.e., the prob that an asthma patient "outgrows" asthma and 
# become undiagnosed with asthma

#iterative method
chosen_province <- "CA"
starting_year <- 1999
end_year <- 2065
stabilization_year <- 2025

get_reassessment_data <- function(
    chosen_province="CA",
    starting_year=1999,
    end_year=2065,
    stabilization_year=2025,
    asthma_max_age=62
){
  
    asthma_inc_model <- read_rds(here("R/asthma_incidence_model.rds"))
    asthma_prev_model <- read_rds(here("R/asthma_prevalence_model.rds"))

  asthma_predictor <- function(age,sex,year,type){
    age <- pmin(age,asthma_max_age)
    
    year <- pmin(year,stabilization_year)
    
    if(type == "prev"){
      return( exp(predict(asthma_prev_model,newdata=data.frame(age,sex,year))) %>% 
                unlist())
    } else{
      return( exp(predict(asthma_inc_model,newdata=data.frame(age,sex,year))) %>% 
                unlist())
    }
    
  }
  
    df_asthma <- expand.grid(age=3:110, sex=c(0,1), year=starting_year:end_year) %>% 
        as.data.frame()
  
    df_asthma <- df_asthma %>% 
        mutate(
            inc=asthma_predictor(age,sex,year,type='inc'),
            prev=asthma_predictor(age,sex,year,type='prev')
        ) %>% 
        # let inc be prev for age = 3
        mutate(inc=ifelse(age==3, prev, inc))
  
    df_asthma_year <- df_asthma %>% group_split(year)
  
    results_assessment <- c()
    years <- c((starting_year+1):end_year)
    ages <- c(4:110)
  
    for( i in 1:(length(years)-1)){
    
        tmp_year <- years[i]
    
        # Get the predicted prevalence for the previous year
        df_prevalence_past <- df_asthma_year[[i]] %>% 
            select(age, year, sex, prev) %>% 
            filter(age!=110) %>%
            pivot_wider(names_from=sex, values_from=prev) %>% 
            select(-age, -year)
    
        # Get the predicted prevalence for the current year
        df_prevalence <- df_asthma_year[[i+1]] %>% 
            select(age, year, sex, prev) %>% 
            filter(age!=3) %>% 
            pivot_wider(names_from=sex, values_from=prev)  %>% 
            select(-age,-year)
    
        # Get the predicted incidence for the current year
        df_incidence <- df_asthma_year[[i+1]] %>% 
            select(age, year, sex, inc) %>% 
            filter(age!=3) %>% 
            pivot_wider(names_from=sex, values_from=inc)  %>% 
            select(-age, -year)

        # estimate tmp_u by assuming P(correct diagnosis) = 1 
        df_assessment <- (df_prevalence - df_incidence * (1 - df_prevalence_past)) / df_prevalence_past
  
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


CA_tuner <- get_reassessment_data(chosen_province="CA",end_year=2066,stabilization_year = 2025)
BC_tuner <- get_reassessment_data(chosen_province="BC",end_year=2043,stabilization_year=2025)

master_assessment <- rbind(CA_tuner,
                           BC_tuner)
master_assessment %>% 
  mutate(M = ifelse(M>1,1,M))
master_assessment <- master_assessment %>% select(year,age,`F`,M,province)

write_csv(master_assessment %>% 
            mutate(M = ifelse(M>1,1,M),
                   `F` = ifelse(`F`>1,1,`F`)),"../src/processed_data/master_asthma_reassessment.csv")


