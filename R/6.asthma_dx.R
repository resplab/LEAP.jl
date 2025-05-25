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
    stabilization_year=2025
){
  
  asthma_inc_model <- read_rds("asthma_incidence_model.rds")
  asthma_prev_model <- read_rds("asthma_prevalence_model.rds")
  asthma_max_age <- 62
  
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
  
  xs <- expand.grid(age=3:110,sex=c(0,1),year=starting_year:end_year) %>% 
    as.data.frame()
  
  df_asthma <- xs %>% 
    mutate(inc = asthma_predictor(age,sex,year,type='inc'),
           prev = asthma_predictor(age,sex,year,type='prev')) %>% 
    # let inc be prev for age = 3
    mutate(inc = ifelse(age==3,prev,inc))
  
  df_asthma_year <- df_asthma %>% 
    group_split(year)
  
  results_assessment <- c()
  years <- c((starting_year+1):end_year)
  ages <- c(4:110)
  
  for( i in 1:(length(years)-1)){
    
    tmp_year <- years[i]
    
    tmp_prev_past <- df_asthma_year[[i]] %>% 
      select(age,year,sex,prev) %>% 
      filter(age!=110) %>%
      pivot_wider(names_from=sex,values_from=prev) %>% 
      select(-age,-year)
    
    tmp_prev <- df_asthma_year[[i+1]] %>% 
      select(age,year,sex,prev) %>% 
      filter(age!=3) %>% 
      pivot_wider(names_from=sex,values_from=prev)  %>% 
      select(-age,-year)
    
    tmp_inc <- df_asthma_year[[i+1]] %>% 
      select(age,year,sex,inc) %>% 
      filter(age!=3) %>% 
      pivot_wider(names_from=sex,values_from=inc)  %>% 
      select(-age,-year)
    
    # prev = prev_past*p(reassessment) + inc*(1-prev_past)
    # p(reassessment) = (prev-inc*(1-prev_past))/prev_past

    # estimate tmp_u by assuming p(correct diagnosis) = 1 
    tmp_assessment <- (tmp_prev - tmp_inc * (1-tmp_prev_past))/tmp_prev_past
    
    # # key step: bound it!
    # tmp_u[tmp_u>1] <- 1
    # tmp_u[tmp_u<0] <- 0
  
    results_assessment[[i]] <- cbind(data.frame(age=4:110,year=tmp_year),tmp_assessment)
  }
  
  results_assessment <- results_assessment %>%
    do.call(rbind,.) %>% 
    as.data.frame()
  
  colnames(results_assessment)[c(3,4)] <- c("F","M")
  
  results_assessment$province <- chosen_province
  
  return(results_assessment)
}


CA_tuner <- get_reassessment_data(chosen_province="CA",end_year=2066,stabilization_year = 2025)
BC_tuner <- get_reassessment_data(chosen_province="BC",end_year=2043,stabilization_year=2025)

# look <- results_assessment %>%
#   do.call(rbind,.) %>%
#   pivot_longer(3:4,names_to="sex",values_to="reassessment") %>%
#   as.data.frame()
# 
# tmp <- gam(log(reassessment)~sex + s(year,age),data=look)
# 
# ggplot(data = look %>% filter(age==12),
#        aes(y=reassessment,x=year)) +
#   geom_line()+
#   ylim(0,1.1)+
#   facet_grid(.~sex)

master_assessment <- rbind(CA_tuner,
                           BC_tuner)
master_assessment %>% 
  mutate(M = ifelse(M>1,1,M))
master_assessment <- master_assessment %>% select(year,age,`F`,M,province)

write_csv(master_assessment %>% 
            mutate(M = ifelse(M>1,1,M),
                   `F` = ifelse(`F`>1,1,`F`)),"../src/processed_data/master_asthma_reassessment.csv")


