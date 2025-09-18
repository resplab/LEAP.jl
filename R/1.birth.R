library(tidyverse)
library(here)

min_year <- 1999

# birth estimate
birth_estimate <- read_csv(here("R/public_dataset","17100005.csv")) %>% 
    filter(GEO %in% c("Canada", "British Columbia")) %>%
    filter(REF_DATE >= min_year) %>%
    filter(`Age group` == "0 years") %>%
    select(1, 2, 4, VALUE) %>%
    rename(
        year=REF_DATE,
        province=GEO,
        sex=Sex,
        N=VALUE
    ) %>%
    mutate(
        province=ifelse(province=="Canada", "CA", "BC"),
        sex=str_sub(sex, 1, 1)
    ) %>% 
    group_by(year, province) %>% 
    mutate(
        max_N=max(N),
        prop=N/max_N
    ) %>%
    ungroup() %>% 
    select(-max_N) %>% 
    filter(sex!='F') %>% 
    select(-prop) %>% 
    pivot_wider(id_cols = 1:2, names_from=sex, values_from=N) %>% 
    mutate(prop_male=M/B) %>% 
    select(-M) %>% 
    select(year,province,B,prop_male) %>% 
    rename(N=B) %>% 
    mutate(projection_scenario='past')


pop_last_year <- max(birth_estimate$year)

pop <- read_csv(here("R/public_dataset", "17100057.csv"))
colnames(pop) <- gsub(' ', '_', colnames(pop))
pop_CA_BC <- pop %>%
    select(REF_DATE, GEO, Projection_scenario, Sex, Age_group, VALUE) %>% 
    mutate(
        Projection_scenario=str_remove(Projection_scenario,"Projection scenario "),
        Projection_scenario=str_remove(Projection_scenario, "\\:.*")) %>% 
    filter(GEO %in% c("Canada", "British Columbia")) %>% 
    filter(REF_DATE>=(pop_last_year + 1)) %>% 
    rename(
        year=REF_DATE,
        province=GEO,
        projection_scenario=Projection_scenario,
        sex=Sex,
        ag=Age_group,
        n=VALUE
    ) %>% 
    mutate(
        province = ifelse(province=="Canada", "CA", "BC"),
        sex = substr(sex,1,1),
        age_remove = str_detect(ag,'to|over|All|Median|Average'),
        age_keep = ag=="100 years and over"
    ) %>% 
    filter(!age_remove | age_keep) %>%
    mutate(
        age=as.numeric(gsub("([0-9]+).*$","\\1",ag)),
        age=ifelse(is.na(age), 0, age)) %>% 
    select(year,sex,age,province,n,projection_scenario) %>%
    mutate(n=n*1000) %>% 
    filter(!is.na(n)) %>% 
    filter(age==0) %>% 
    pivot_wider(names_from=sex, values_from=n) %>% 
    mutate(prop_male=M/B) %>% 
    select(year, province, B, prop_male, projection_scenario) %>% 
    rename(N=B)

birth_estimate <- rbind(birth_estimate,pop_CA_BC)

write_csv(birth_estimate, here("src/processed_data/master_birth_estimate.csv"))
