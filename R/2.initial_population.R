library(tidyverse)
library(here)

min_year <- 2000

# initial population ------------------------------------------------------
initial_pop <- read_csv(here("R/public_dataset/17100005.csv")) %>% 
    filter(GEO %in% c('Canada', "British Columbia")) %>% 
    filter(REF_DATE >= min_year) %>% 
    mutate(
        sex=str_sub(Sex,1,1),
        age_remove = str_detect(`Age group`,'to|over|All|Median|Average'),
        age_keep = `Age group`=="100 years and over"
    )

initial_pop_distribution <- initial_pop %>% 
    filter(!age_remove | age_keep) %>%
    mutate(age =as.numeric(gsub("([0-9]+).*$","\\1",`Age group`))) %>% 
    rename(n=VALUE) %>% 
    select(1,sex,age,GEO,n) %>% 
    rename(
        year=REF_DATE,
        province=GEO
    ) %>% 
    mutate(
        province=case_when(
            province=="Canada" ~ "CA",
            province=="British Columbia" ~ "BC",
            TRUE ~ "NONE")
    ) %>% 
    filter(sex!="B")

# approximate
tmp <- initial_pop_distribution %>% 
    filter(is.na(n))

tmp2 <- initial_pop_distribution %>% 
    filter(year %in% ((tmp$year %>% unique())+1)) %>% 
    filter(age %in% ((tmp$age %>% unique()) +1)) %>% 
    mutate(age=age-1)

tmp3 <- tmp %>% 
    select(-n) %>% 
    left_join(
        tmp2 %>% select(-year),
        by=c('sex','age','province')
    ) %>% 
    select(n) %>% 
    unlist()

tmp3 <- initial_pop_distribution$n[which(is.na(initial_pop_distribution$n))]

initial_pop_distribution <- initial_pop_distribution %>% 
    filter(!is.na(n)) %>% 
    group_by(year,province) %>% 
    mutate(total_n = sum(n)) %>% 
    mutate(prop=n/total_n)

initial_pop_distribution_both <- initial_pop_distribution %>% 
    group_by(year,age,province) %>% 
    mutate(n_both=sum(n)) %>% 
    mutate(prop_male = n/n_both) %>% 
    ungroup() %>% 
    filter(sex=='M') %>% 
    select(year, age, province,prop_male)

initial_pop_distribution_prop <- initial_pop_distribution %>% 
    group_by(year,age,province) %>% 
    summarise(n=sum(n)) %>% 
    ungroup() %>% 
    arrange(year,province)

initial_pop_distribution_birth_n <- initial_pop_distribution_prop %>% 
    filter(age==0) %>% 
    select(-age) %>% 
    rename(n_birth=n)

initial_pop_distribution_prop_past <- initial_pop_distribution_prop %>% 
    left_join(initial_pop_distribution_birth_n, by=c("year","province")) %>% 
    mutate(prop=n/n_birth) %>% 
    left_join(initial_pop_distribution_both, by=c("year","age","province"))


initial_pop <- read_csv(here("R/public_dataset","17100057.csv"))
colnames(initial_pop) <- gsub(' ', '_',colnames(initial_pop))
pop_projection_BC_CA <- initial_pop %>% 
    select(REF_DATE,GEO,Projection_scenario,Sex,Age_group,VALUE) %>% 
    mutate(
        Projection_scenario=str_remove(Projection_scenario, "Projection scenario "),
        Projection_scenario=str_remove(Projection_scenario, "\\:.*")
    ) %>% 
    filter(GEO %in% c('Canada', "British Columbia"))

write_rds(pop_projection_BC_CA, here("src/processed_data/pop_projection_BC_CA.rds"))


initial_pop_distribution <- pop_projection_BC_CA %>%
    filter(REF_DATE >= max(initial_pop_distribution_prop$year)) %>%
    mutate(
        Projection_scenario = str_remove(Projection_scenario,"Projection scenario "),
        Projection_scenario = str_remove(Projection_scenario, "\\:.*")
    ) %>%
    mutate(sex=str_sub(Sex,1,1),
            age_remove = str_detect(Age_group,'to|over|All|Median|Average'),
            age_keep = Age_group=="100 years and over") %>% 
    filter(!age_remove | age_keep) %>%
    mutate(age = as.numeric(gsub("([0-9]+).*$","\\1",Age_group)),
            age = ifelse(is.na(age),0,age)) %>% 
    rename(n=VALUE,
            projection_scenario = Projection_scenario ) %>% 
    select(1,projection_scenario,sex,age,GEO,n) %>% 
    rename(year=REF_DATE,
            province=GEO) %>% 
    mutate(n=n*1000) %>% 
    mutate(province = case_when(province=="Canada" ~ "CA",
                                province=="British Columbia" ~ "BC",
                                TRUE ~ "NONE")) %>% 
    filter(sex!="B") %>% 
    filter(!is.na(n))



# initial_pop <- readRDS(here("public_dataset","pop_projection_BC_CA.rds")) %>% 
#   filter(GEO %in% c('Canada',"British Columbia")) %>% 
#   filter(REF_DATE >= max(initial_pop_distribution_prop$year)) %>% 
#   mutate(Projection_scenario = str_remove(Projection_scenario,"Projection scenario "),
#          Projection_scenario = str_remove(Projection_scenario, "\\:.*")) %>% 
#   mutate(sex=str_sub(Sex,1,1),
#          age_remove = str_detect(Age_group,'to|over|All|Median|Average'),
#          age_keep = Age_group=="100 years and over")
# 
# initial_pop_distribution <- initial_pop %>% 
#   filter(!age_remove | age_keep) %>%
#   mutate(age = as.numeric(gsub("([0-9]+).*$","\\1",Age_group)),
#          age = ifelse(is.na(age),0,age)) %>% 
#   rename(n=VALUE,
#          projection_scenario = Projection_scenario ) %>% 
#   select(1,projection_scenario,sex,age,GEO,n) %>% 
#   rename(year=REF_DATE,
#          province=GEO) %>% 
#   mutate(n=n*1000) %>% 
#   mutate(province = case_when(province=="Canada" ~ "CA",
#                               province=="British Columbia" ~ "BC",
#                               TRUE ~ "NONE")) %>% 
#   filter(sex!="B") %>% 
#   filter(!is.na(n))

# note that for province, projections are not available 2044 onwards

initial_pop_distribution %>% 
  group_by(year,province,projection_scenario) %>% 
  mutate(total_n = sum(n)) %>% 
  mutate(prop=n/total_n) -> initial_pop_distribution

initial_pop_distribution_both <- initial_pop_distribution %>% 
    group_by(year,age,province,projection_scenario) %>% 
    mutate(n_both=sum(n)) %>% 
    mutate(prop_male = n/n_both) %>% 
    ungroup() %>% 
    filter(sex=='M') %>% 
    select(year,age,province,prop_male,projection_scenario)

initial_pop_distribution_prop <- initial_pop_distribution %>% 
    group_by(year,age,province,projection_scenario) %>% 
    summarise(n=sum(n)) %>% 
    ungroup() %>% 
    arrange(year,province,projection_scenario)

initial_pop_distribution_birth_n <- initial_pop_distribution_prop %>% 
    filter(age==0) %>% 
    select(-age) %>% 
    rename(n_birth=n)

initial_pop_distribution_prop_project <- initial_pop_distribution_prop %>% 
    left_join(initial_pop_distribution_birth_n, by=c("year","province", "projection_scenario")) %>% 
    mutate(prop=n/n_birth) %>% 
    left_join(initial_pop_distribution_both, by=c("year", "age", "province", "projection_scenario")) 

master_initial_pop_distribution_prop <- initial_pop_distribution_prop_past %>% 
    mutate(projection_scenario="past") %>% 
    rbind(.,initial_pop_distribution_prop_project %>% 
    select(year,age,province,n,n_birth,prop,prop_male,projection_scenario))

write_csv(
    master_initial_pop_distribution_prop,
    here('src/processed_data/master_initial_pop_distribution_prop.csv')
)