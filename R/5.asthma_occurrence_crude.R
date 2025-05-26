# source: https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1310009608

library(tidyverse)
library(here)
library(mgcv)
max_year <- 2019
df_raw <- read_csv(here("public_dataset","13100096.csv"))
chosen_province <- "British Columbia"
baseline_year <- 2000
# chosen_province <- "Canada"


load_asthma_df_admin <- function(starting_year=baseline_year) {
    df <- readxl::read_xlsx(here("R/private_dataset/asthma_inc_prev.xlsx"), sheet=1)

    df <- df %>% filter(age_group_desc != "<1 year") %>% 
        mutate(year=substr(fiscal_year,1,4) %>% as.numeric()) %>% 
        rename(sex=gender) %>% 
        filter(year>=starting_year) %>% 
        rename(age_group = age_group_desc)

    lapply(df$age_group,function(x){
        ceiling(mean(parse_number(str_split(x,"-")[[1]])))
    })  %>% unlist() -> df$age

    # Key assumption: set the incidence level at age 3 to the prevelance level
    df <- df %>% 
        mutate(incidence = ifelse(age==3, prevalence, incidence)) %>% 
        mutate(age= ifelse(age==90, 100, age))

    return(df)
}

df_admin <- load_asthma_df_admin()

# BC asthma prev and inc equation -----------------------------------------

load_asthma_df_bc <- function(starting_year=baseline_year) {
    df <- readxl::read_xlsx(here("R/private_dataset/asthma_inc_prev.xlsx"), sheet=1)

    df <- df %>% filter(age_group_desc != "<1 year") %>% 
        mutate(year=substr(fiscal_year, 1, 4) %>% as.numeric()) %>% 
        rename(sex=gender) %>% 
        filter(year>=starting_year) %>% 
        rename(age_group = age_group_desc) %>% 
        mutate(prev_sd = sqrt(prevalence_numerator)*qnorm(0.975),
                prev_upper = (prevalence_numerator+prev_sd)/pop,
                prev_lower = (prevalence_numerator-prev_sd)/pop)

    lapply(df$age_group,function(x){
        ceiling(mean(parse_number(str_split(x,"-")[[1]])))
    })  %>% 
    unlist() -> df$age

    # Key assumption: asthma starts at age 3
    # Set the inc = prev at age 3
    df <- df %>% 
        mutate(
            incidence=ifelse(age==3,prevalence,incidence),
            incidence_numerator=ifelse(age==3,prevalence_numerator,incidence_numerator)
        )

    return(df)
}

master_BC_asthma <- load_asthma_df_bc()


# INCIDENCE MODEL: BC ---------------------------------------------------------
# INCIDENCE: log(incidence) ~ sex(year + sex*poly(age,5))

generate_incidence_model <- function(df_asthma, min_age=3, max_age=65, starting_year=baseline_year) {
    # Create incidence dataframe
    df <- df_asthma %>% select(year, sex, age, incidence)

    # Filter for year <= max_year
    df <- df %>% filter(year <= max_year)
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

    chosen_year <- seq(2000,2020,by=2)

    ggplot(data=df_inc %>% 
            filter(year %in% chosen_year) %>% 
            mutate(year=as.factor(year)),
        aes(x=age,y=y,color=year)) +
    geom_line() +
    geom_line(data=inc_pred %>% 
                filter(year %in% chosen_year) %>% 
                mutate(year=as.factor(year)),aes(x=age,y=y,color=year),
                linetype="dashed") +
    ylab("Asthma incidence per 100 in BC") +
    facet_grid(.~sex)

    write_rds(model, here("R/private_dataset/asthma_incidence_model.rds"))
}

generate_incidence_model(df_asthma=master_BC_asthma, min_age=3, max_age=65)

# PREV: sex*poly(year,degree=2)*poly(age,degree=5)

df_prev <- master_BC_asthma %>% 
  select(year,sex,age,prevalence) %>% 
  filter(year >= baseline_year & year <= max_year) %>%
  mutate(sex=as.numeric(sex=="M"),
         y= prevalence) %>% 
  filter(age <= 65)


prev_model <- gam(log(prevalence)~ sex*poly(year,degree=2)*poly(age,degree=5),
                  data=df_prev)

summary(prev_model)

prev_pred <- expand.grid(year=c(2000:2065),sex=c(0,1),age=seq(3,62,by=1)) 

prev_pred$y <- exp(predict(prev_model,newdata=prev_pred))

df_prev_pred <- df_prev %>% 
  left_join(prev_pred,by=c("year","sex","age"))

chosen_year <- seq(2000,2025,by=2)

ggplot(data=df_prev %>% 
         filter(year %in% chosen_year) %>% 
         mutate(year=as.factor(year)),
       aes(x=age,y=y,color=year)) +
  geom_line() +
  geom_line(data=prev_pred %>% 
              filter(year %in% chosen_year) %>% 
              mutate(year=as.factor(year)),aes(x=age,y=y,color=year),
            linetype="dashed") +
  ylab("Asthma prevalence per 100 in BC") +
  facet_grid(.~sex)

# test poly

z <- poly(df_prev$year,2)

z_attr <- attributes(z)
alpha <- z_attr$coefs$alpha
nd <- z_attr$coefs$norm2
degree <- 5

basis <- function(x,alpha,nd,degree){
  fs <-  1/sqrt(nd[2])
  fs <- append(fs,(x-alpha[1]) / sqrt(nd[3]))
  if(degree>1){
    for(i in 2:degree){
      fs <- append(fs,((x-alpha[i]) * sqrt(nd[i+1]) * fs[i] - nd[i+1] / sqrt(nd[i]) * fs[i-1]) / sqrt(nd[i+2]))
    }
  }
  fs[-1]
}

write_rds(prev_model,"asthma_prevalence_model.rds")




# generate crude asthma inc & prev  ----------------------------------------

library(mgcv)
library(tidyverse)
prev_model <- read_rds("asthma_prevalence_model.rds")
inc_model <- read_rds("asthma_incidence_model.rds")
stabilization_year <- 2025
max_age <- 63
df <- expand.grid(year=2000:2065,sex=c(0:1),age=3:110) %>% 
  as.data.frame() %>% 
  mutate(prev = exp(predict(prev_model,data.frame(year=pmin(2025,year),sex,age=pmin(age,max_age)))),
         inc = exp(predict(inc_model,data.frame(year=pmin(2025,year),sex,age=pmin(age,max_age))))) %>% 
  mutate(prev=as.numeric(prev),
         inc = as.numeric(inc))

ggplot(data=df %>% 
         mutate(sex = ifelse(sex==1,"Male","Female")) %>% 
         filter(year %in% seq(2000,2025,5)) %>% 
         mutate(year=as.factor(year)),aes(x=age,y=prev,col=year)) +
  geom_line() +
  xlim(c(0,60)) +
  facet_grid(.~sex) +
  ylab("Crude asthma prevalence (per 100)")+
  xlab("Age (year)")+
  theme_classic(base_size=20) +
  theme(legend.position='top',
        legend.title=element_blank())

ggplot(data=df %>% 
         mutate(sex = ifelse(sex==1,"Male","Female")) %>% 
         filter(year %in% seq(2000,2025,5)) %>% 
         mutate(year=as.factor(year)),aes(x=age,y=inc,col=year)) +
  geom_line() +
  xlim(c(0,60)) +
  facet_grid(.~sex) +
  ylab("Crude asthma incidence (per 100)")+
  xlab("Age (year)")+
  theme_classic(base_size=20) +
  theme(legend.position='top',
        legend.title=element_blank())

# write_csv(df,"master_asthma_prev_inc.csv")

