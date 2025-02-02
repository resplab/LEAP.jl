library(tidyverse)
library(parallel)

baseline_year <- 2000
last_year <- 2020
projected_last_year <- 2068
life_table_list <- c()
immigration_list <- c()
emigration_list <- c()
life_table_checker <- c()

provinces <- c("BC","CA")
desired_life_expectancys <- list(c(84.6,88.0),c(87,90.1)) # BC: male, female; CANADA: male, female
calibration_years <- c(2043, 2068)

for(province_index in 1:length(provinces)) {
  chosen_province <- provinces[province_index]
  calibration_year <- calibration_years[province_index]
  desired_life_expectancy <- desired_life_expectancys[[province_index]]
  print(chosen_province)
  
# life table --------------------------------------------------------------

life_table <- read_csv("life_table.csv") %>% 
  filter(province==chosen_province)

death_final_year <- max(life_table$year)

ref_life_table <- life_table %>% 
  filter(year==death_final_year)

death_adjustment <- function(p,year,beta){
        p <- pmin(p, 0.9999999999)
        odds <- p/(1 - p)*exp(year*beta)
  pmax(pmin(odds/(1+odds),1),0)
}

    project_life_year <- function(
        ref_lt, proj_year, baseyear, beta_year
    ){
        tmp <- death_adjustment(ref_lt$prob_death,  proj_year-baseyear,beta_year)
  tmp_factor <- tmp/ref_lt$prob_death
  ref_lt %>% 
    mutate(prob_death = tmp_factor*prob_death,
           se = tmp_factor*se,
           year = proj_year)
}

projected_life_table <- c()

    obj_function_female <- function(mortality_year_adjustment, SEX="F"){

  for(yr in 1:(projected_last_year-death_final_year)){
            projected_life_table[[yr]] <- project_life_year(
                ref_life_table, yr+death_final_year, baseyear=death_final_year,
                beta_year = mortality_year_adjustment
            )
  }
  
  life_expectancy_calculator <- function(lf){
    lf$I <- NA
    lf$I[1] <- 100000
      for(i in 2:nrow(lf)){
        lf$I[i] <- lf$I[i-1]*(1-lf$q[i-1])
      }
      
            lf %>% mutate(d = I*q, L = lead(I)+0.5*d) -> lf
    lf$L[1] <- lf$L[2] + 0.1*lf$d[1]
    lf$L[111] <- lf$I[111]*1.4
      
    lf$TT <- rev(cumsum(rev(lf$L)))
            lf %>% mutate(E = TT/I) -> sol
      sol$E[1]
  }
  
  
  projected_life_table <- do.call(rbind,projected_life_table)
  
  projected_life_table %>% 
    filter(sex==SEX & year==calibration_year) %>% 
    select(age,prob_death) %>% 
    rename(q = prob_death) -> lf
  
  return(life_expectancy_calculator(lf) - desired_life_expectancy[as.numeric(SEX=="F")+1])
}


obj_function_male <- function(mortality_year_adjustment,SEX="M"){
  
  for(yr in 1:(projected_last_year-death_final_year)){
    projected_life_table[[yr]] <- project_life_year(ref_life_table,
                                                    yr+death_final_year,
                                                    baseyear=death_final_year,
                                                    beta_year = mortality_year_adjustment)
  }
  
  life_expectancy_calculator <- function(lf){
    lf$I <- NA
    lf$I[1] <- 100000
    for(i in 2:nrow(lf)){
      lf$I[i] <- lf$I[i-1]*(1-lf$q[i-1])
    }
    
    lf %>% 
            mutate(d = I*q, L = lead(I)+0.5*d) -> lf
    lf$L[1] <- lf$L[2] + 0.1*lf$d[1]
    lf$L[111] <- lf$I[111]*1.4
    
    lf$TT <- rev(cumsum(rev(lf$L)))
            lf %>% mutate(E = TT/I) -> sol
    sol$E[1]
  }
  
  projected_life_table <- do.call(rbind,projected_life_table)
  
  projected_life_table %>% 
    filter(sex==SEX & year==calibration_year) %>% 
    select(age,prob_death) %>% 
    rename(q = prob_death) -> lf
  
  return(life_expectancy_calculator(lf) - desired_life_expectancy[as.numeric(SEX=="F")+1])
}

projected_life_table_male <- projected_life_table_female <- c()

for(yr in 1:(projected_last_year-death_final_year)){
        projected_life_table_male[[yr]] <- project_life_year(
            ref_life_table,
                                                  yr+death_final_year,
                                                  baseyear=death_final_year,
            beta_year=  uniroot(obj_function_male, interval = c(-0.03,-0.01),tol = 0.00001)$root
        )
} 

projected_life_table_male <- do.call(rbind,projected_life_table_male) %>% 
  filter(sex=="M")

for(yr in 1:(projected_last_year-death_final_year)){
        projected_life_table_female[[yr]] <- project_life_year(
            ref_life_table,
                                                       yr+death_final_year,
                                                       baseyear=death_final_year,
            beta_year=uniroot(obj_function_female,interval = c(-0.03,-0.01),tol = 0.00001)$root
        )

} 

projected_life_table_female <- do.call(rbind,projected_life_table_female) %>% 
  filter(sex=="F")

projected_life_table <- rbind(projected_life_table_male,projected_life_table_female)

life_table <- rbind(life_table,projected_life_table)
life_table_list <- rbind(life_table_list,life_table)

# pop growth --------------------------------------------------------------

pop_est <- read_csv("master_initial_pop_distribution_prop.csv") %>% 
  filter(province==chosen_province) %>%
  filter(year >= baseline_year) %>% 
  mutate(Male = prop_male*n,
         Female = (1-prop_male)*n) %>% 
  select(year,age,province,Male,Female,projection_scenario) %>% 
  pivot_longer(4:5,names_to="sex",values_to="n") %>% 
  mutate(sex=ifelse(sex=="Male","M","F")) %>% 
  select(year,sex,age,province,n,projection_scenario) %>% 
  arrange(year,desc(sex),age,province,projection_scenario)

pop_scenarios <- pop_est$projection_scenario %>% unique()
pop_scenarios <- pop_scenarios[-which(pop_scenarios=="past")]

    max_pop_year <- min(max(pop_est$year), 2065)

for(i in 1:length(pop_scenarios)){
  tmp_pop <- pop_est %>% 
    filter(projection_scenario %in% c("past",pop_scenarios[i])) %>% 
    filter(!(year==2021 & projection_scenario %in% c(pop_scenarios[i]))) %>% 
    select(-projection_scenario)
  
  df_diff <- expand.grid(year=(baseline_year+1):max_pop_year,
                         age = 1:100,
                         sex = c("F","M")) %>% 
    mutate(n = 0)
  

  tmp_combined <- tmp_pop %>% 
    left_join(life_table,by=c("age",'sex','year','province'))
  
  mclapply(X=1:nrow(df_diff),mc.cores = 7,FUN = function(j){
    YEAR <- df_diff$year[j]
    AGE <- df_diff$age[j]
    SEX <- df_diff$sex[j]
    tmp <- tmp_combined %>% 
    filter((year %in% c(YEAR,YEAR-1) & age %in% c(AGE,AGE-1) & sex==SEX) )
    tmp$n[nrow(tmp)] - tmp$n[1]*(1-tmp$prob_death[1])
  } ) -> tmp_n
  
  df_diff$n <- unlist(tmp_n)
  
  df_diff %>% 
    arrange(year,age,sex) %>%
    filter(age<=100) -> look
  
  tmp_pop_birth <- tmp_pop %>% 
    filter(age==0) %>% 
    select(1,2,5) %>% 
    rename(n_birth = n) %>% 
    group_by(year) %>% 
    summarise(n_birth=sum(n_birth))
  
  immigration_n <- look %>% 
    mutate(n= ifelse(n<0,0,n)) %>% 
    left_join(tmp_pop_birth,by=c("year")) %>% 
    mutate(n_prop_birth = n/n_birth) %>% 
    select(year,age,sex,n_prop_birth) %>% 
    group_by(year) %>% 
    mutate(tot=sum(n_prop_birth)) %>% 
    ungroup() %>% 
    mutate(weights = n_prop_birth/tot) %>% 
    select(-tot) %>% 
    mutate(sex=as.numeric(sex=="M")) %>% 
            mutate(province=chosen_province, proj_scenario = pop_scenarios[i])
  
  immigration_list <- rbind(immigration_list,immigration_n)
  
  emigration_n <- look %>% 
    mutate(n= ifelse(n>0,0,-n)) %>% 
    left_join(tmp_pop,by=c("year",'sex','age')) %>% 
    mutate(prob = n.x/n.y) %>% 
    select(year,age,sex,n.x,prob) %>% 
    rename(n=n.x) %>% 
    select(year,age,sex,prob) %>% 
    pivot_wider(names_from=sex,values_from=prob) %>%
            mutate(province=chosen_province, proj_scenario = pop_scenarios[i])
  
  emigration_list <- rbind(emigration_list,emigration_n)
    }
}

write_csv(life_table_list,here("../src","processed_data","master_life_table.csv"))
write_csv(immigration_list,here("../src","processed_data","master_immigration_table.csv"))
write_csv(emigration_list,here("../src","processed_data","master_emigration_table.csv"))

# offset
immigration_list <- read_csv("../src/processed_data/master_immigration_table.csv")

immigration_list %>% 
  filter(year==2023) %>% 
  filter(proj_scenario=="M3") %>% 
  filter(province=="CA") -> tmp

tmp$n_prop_birth[1:2] <- tmp$n_prop_birth[1:2]*4
tmp$weights <- tmp$n_prop_birth/sum(tmp$n_prop_birth)

immigration_list %>% 
  filter(!(year==2023 & proj_scenario=="M3" & province=="CA")) -> look

rbind(look,tmp) %>% 
  arrange(year,age,sex,province,proj_scenario) -> final_look

write_csv(final_look,"../src/processed_data/master_immigration_table_modified.csv")
