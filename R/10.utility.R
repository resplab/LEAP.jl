library(tidyverse)
library(janitor)
library(here)

eq5d <- readxl::read_xlsx("public_dataset/canada_eq_5d.xlsx", 1) %>%
    janitor::clean_names() %>%
    pivot_longer(-c(1:3), values_to="eq5d", names_to="sex") %>%
    mutate(sex=ifelse(sex == "male", 1, 0)) %>%
    rowwise() %>%
    mutate(
        se=str_extract(eq5d, "(?<=\\().*(?=\\))") %>% as.numeric()
    ) %>%
    mutate(
        eq5d=str_replace(eq5d, " \\s*\\([^\\)]+\\)", "") %>% as.numeric()
    )

results <- c()

counter <- 1
for (tmp_age in 18:111){
    for (tmp_sex in c(0, 1)){
        tmp_result <- eq5d %>%
            filter(tmp_age <= age_upper & tmp_age >= age_lower) %>%
            filter(sex == tmp_sex)
        results[[counter]] <- data.frame(
            age=tmp_age, sex=tmp_sex, eq5d=tmp_result$eq5d, se=tmp_result$se
        )
        counter <- counter + 1
    }
}

results <- results %>%
    do.call(rbind, .) %>%
    as.data.frame()

min_age <- results$age %>% min()
tmp_age <- seq(0, min_age - 1, by=1)

# linearly interpolate btw ages 0 and 18
for (tmp_sex in c(0, 1)){
    slope <- (
        1 - (
            results %>% filter(sex == tmp_sex & age == min_age) %>% select(eq5d) %>% unlist()
        )
    ) / (min_age)
    results <- rbind(
        results,
        data.frame(age=tmp_age, sex=tmp_sex, eq5d=1 - slope * tmp_age, se=0)
    )
}

eq5d_canada <- results %>%
    arrange(sex, age)

write_csv(eq5d_canada, here("src/processed_data/eq5d_canada.csv"))
