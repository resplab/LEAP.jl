library(tidyverse)
library(parallel)
library(here)

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

get_prob_death_projected = function(prob_death, year_index, beta) {
    prob_death = pmin(prob_death, 0.9999999999)
    odds = (prob_death/(1 - prob_death)) * exp(year_index * beta)
    return(pmax(pmin(odds/(1 + odds), 1), 0))
}


get_projected_life_table_single_year <- function(
    ref_life_table, death_final_year, year_index, beta_year
){
    prob_death_projected = get_prob_death_projected(
        prob_death=ref_life_table$prob_death,
        year_index=year_index,
        beta=beta_year
    )
    return(
        ref_life_table %>% 
            mutate(
                prob_death=prob_death_projected,
                se=se * prob_death_projected / ref_life_table$prob_death,
                year=year_index + death_final_year
            )
    )
}


life_expectancy_calculator <- function(life_table_year){
    life_table_year$I <- NA
    life_table_year$I[1] <- 100000
    for(i in 2:nrow(life_table_year)){
        life_table_year$I[i] <- life_table_year$I[i-1]*(1-life_table_year$q[i-1])
    }
    
    life_table_year <- life_table_year %>% mutate(d=I*q, L=lead(I) + 0.5*d)
    life_table_year$L[1] <- life_table_year$L[2] + 0.1*life_table_year$d[1]
    life_table_year$L[111] <- life_table_year$I[111]*1.4
    
    life_table_year$T <- rev(cumsum(rev(life_table_year$L)))
    life_table_year = life_table_year %>% mutate(E=T/I)
    return(life_table_year$E[1])
}


beta_year_optimizer <- function(
    mortality_year_adjustment, SEX, projected_last_year, death_final_year,
    projected_life_table, ref_life_table, desired_life_expectancy, calibration_year
){

    for(yr in 1:(projected_last_year - death_final_year)){
        projected_life_table[[yr]] <- get_projected_life_table_single_year(
            ref_life_table=ref_life_table,
            death_final_year=death_final_year,
            year_index=yr,
            beta_year=mortality_year_adjustment
        )
    }

    projected_life_table <- do.call(rbind, projected_life_table)
    
    lf <- projected_life_table %>% 
        filter(sex==SEX & year==calibration_year) %>% 
        select(age, prob_death) %>% 
        rename(q=prob_death)

    life_expectancy = life_expectancy_calculator(lf)
    return(life_expectancy - desired_life_expectancy[as.numeric(SEX=="F") + 1])
}


get_prev_year_population <- function(row, tmp_combined){
    YEAR <- row$year
    AGE <- row$age
    SEX <- row$sex
    tmp <- tmp_combined %>% 
        filter((year %in% c(YEAR, YEAR-1) & age %in% c(AGE, AGE-1) & sex==SEX))
    tmp$n[nrow(tmp)] - tmp$n[1]*(1-tmp$prob_death[1])
}

for(province_index in 1:length(provinces)) {
    chosen_province <- provinces[province_index]
    calibration_year <- calibration_years[province_index]
    desired_life_expectancy <- desired_life_expectancys[[province_index]]
    print(chosen_province)
  
    # life table --------------------------------------------------------------

    life_table <- read_csv(here("life_table.csv")) %>% 
    filter(province==chosen_province)

    death_final_year <- max(life_table$year)

    ref_life_table <- life_table %>% filter(year==death_final_year)



    projected_life_table <- c()

    projected_life_table_male <- projected_life_table_female <- c()

    for(yr in 1:(projected_last_year-death_final_year)){
        beta_year <- uniroot(
            beta_year_optimizer,
            SEX="M",
            projected_last_year=projected_last_year,
            death_final_year=death_final_year,
            projected_life_table=projected_life_table,
            ref_life_table=ref_life_table,
            desired_life_expectancy=desired_life_expectancy,
            calibration_year=calibration_year,
            interval=c(-0.03,-0.01),
            tol=0.00001
        )$root
        projected_life_table_male[[yr]] <- get_projected_life_table_single_year(
            ref_life_table=ref_life_table,
            death_final_year=death_final_year,
            year_index=yr,
            beta_year=beta_year
        )
    } 

    projected_life_table_male <- do.call(rbind, projected_life_table_male) %>% filter(sex=="M")

    for(yr in 1:(projected_last_year-death_final_year)){
        beta_year <- uniroot(
            beta_year_optimizer,
            SEX="F",
            projected_last_year=projected_last_year,
            death_final_year=death_final_year,
            projected_life_table=projected_life_table,
            ref_life_table=ref_life_table,
            desired_life_expectancy=desired_life_expectancy,
            calibration_year=calibration_year,
            interval=c(-0.03,-0.01),
            tol=0.00001
        )$root
        projected_life_table_female[[yr]] <- get_projected_life_table_single_year(
            ref_life_table=ref_life_table,
            death_final_year=death_final_year,
            year_index=yr,
            beta_year=beta_year
        )
    } 

    projected_life_table_female <- do.call(rbind, projected_life_table_female) %>% filter(sex=="F")

    projected_life_table <- rbind(projected_life_table_male, projected_life_table_female)

    life_table <- rbind(life_table, projected_life_table)
    life_table_list <- rbind(life_table_list, life_table)

    # pop growth --------------------------------------------------------------

    df_population <- read_csv(here("src/processed_data", "master_initial_pop_distribution_prop.csv"))

    # Select province and years
    df_population <- df_population %>% 
        filter(province==chosen_province) %>%
        filter(year >= baseline_year)

    # Get the total number of male / female population for given year/age/projection_scenario
    df_population <- df_population %>% 
        mutate(M=prop_male*n, F=(1-prop_male)*n)
        
    df_population <- df_population %>% 
        select(year, age, province, M, F, projection_scenario)
        
    df_population <- df_population %>% 
        pivot_longer(4:5, names_to="sex", values_to="n")
        
    df_population <- df_population %>%
        select(year, sex, age, province, n, projection_scenario) %>% 
        arrange(year, desc(sex), age, province, projection_scenario)

    pop_scenarios <- df_population$projection_scenario %>% unique()
    pop_scenarios <- pop_scenarios[-which(pop_scenarios=="past")]

    max_pop_year <- min(max(df_population$year), 2065)

    for(i in 1:length(pop_scenarios)){
        df_proj <- df_population %>% 
            filter(projection_scenario %in% c("past", pop_scenarios[i])) %>% 
            filter(!(year==2021 & projection_scenario %in% c(pop_scenarios[i]))) %>% 
            select(-projection_scenario)
  
        df_diff <- expand.grid(
            year=(baseline_year+1):max_pop_year,
            age=1:100,
            sex=c("F","M")
        ) %>% mutate(n = 0)
  
        df_proj <- df_proj %>% left_join(life_table, by=c("age",'sex','year','province'))
  
        delta_n <- mclapply(
            X=split(df_diff, 1:nrow(df_diff)),
            mc.cores=7,
            FUN=get_prev_year_population,
            tmp_combined=df_proj
        )
  
        df_diff$delta_n <- unlist(delta_n)
  
        df_diff <- df_diff %>%
            arrange(year, age, sex) %>%
            filter(age<=100)
  
        df_birth <- df_proj %>% 
            filter(age==0) %>% 
            select(1, 2, 5) %>% 
            rename(n_birth=n) %>% 
            group_by(year) %>% 
            summarise(n_birth=sum(n_birth))
  
        # If delta_n is negative, set n = 0, else n = delta_n
        df_immigration <- df_diff %>% 
            mutate(n_immigrants=ifelse(delta_n < 0, 0, delta_n))

        # Add the n_birth column
        df_immigration <- df_immigration %>%
            left_join(df_birth, by=c("year"))
            
        # Get the proportion of immigrants relative to the number of people born that year
        df_immigration <- df_immigration %>% 
            mutate(prop_immigrants_birth = n_immigrants / n_birth)

        df_immigration <- df_immigration %>%
            select(year, age, sex, prop_immigrants_birth)
        
        df_immigration <- df_immigration %>%
            group_by(year) %>% 
            mutate(tot=sum(prop_immigrants_birth)) %>% 
            ungroup() %>% 
            mutate(weights = prop_immigrants_birth/tot) %>% 
            select(-tot) %>% 
            mutate(province=chosen_province, proj_scenario = pop_scenarios[i])
  
        immigration_list <- rbind(immigration_list, df_immigration)

        df_emigration <- df_diff %>% 
            mutate(n= ifelse(n>0,0,-n)) %>% 
            left_join(df_proj,by=c("year",'sex','age')) %>% 
            mutate(prob = n.x/n.y) %>% 
            select(year,age,sex,n.x,prob) %>% 
            rename(n=n.x) %>% 
            select(year,age,sex,prob) %>% 
            pivot_wider(names_from=sex,values_from=prob) %>%
            mutate(province=chosen_province, proj_scenario = pop_scenarios[i])

        emigration_list <- rbind(emigration_list, df_emigration)
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
