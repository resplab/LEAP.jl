library(tidyverse)
library(here)
library(mgcv)
library(roptim)
source(here("R/calibration_helper_function.R"))
options(dplyr.summarise.inform = FALSE)

chosen_province <- "CA"
max_cal_year <- 2065 # 2065 for CA; 2043 for BC
min_cal_year <- 2000
stabilization_year <- 2025
MAX_ASTHMA_AGE <- 62
MIN_ASTHMA_AGE <- 3
# odds ratio between asthma prevalence at age 3 and family history (CHILD Study)
OR_ASTHMA_AGE_3 <- 1.13
# odds ratio between asthma prevalence at age 5 and family history (CHILD Study)
OR_ASTHMA_AGE_5 <- 2.4
INC_BETA_PARAMS <- c((log(OR_ASTHMA_AGE_5) - log(OR_ASTHMA_AGE_3)) / 2, -0.225)
# the probability that one or more parents have asthma (CHILD Study)
PROB_FAM_HIST <- 0.2927242


asthma_predictor <- function(age, sex, year, type, asthma_inc_model, asthma_prev_model) {

    age <- pmin(age, MAX_ASTHMA_AGE)
    year <- pmin(year, stabilization_year)
  
    if(type=="prev"){
        return(exp(predict(asthma_prev_model, newdata=data.frame(age, sex, year))) %>% 
                unlist())
    } else {
        return(exp(predict(asthma_inc_model, newdata=data.frame(age, sex, year))) %>% 
                unlist())
    }
}

#' Compute the probability of number of courses of antibiotics during infancy.
#' 
#' @param chosen_year The birth year of the infant.
#' @param chosen_sex The sex of the infant; 0 = female, 1 = male.
#' @param model_abx The fitted Negative Binomial model for the number of courses of antibiotics.
#' @returns A dataframe with the probability of the number of courses of antibiotics,
#' ranging from 0 - 5+.
p_antibiotic_exposure <- function(chosen_year, chosen_sex, model_abx){
    # 2025 for females
    # 2028 for males
    # to cap it
    if(chosen_sex == 1){
        chosen_year <- min(2028 - 1, chosen_year)
    }  else{
        chosen_year <- min(2025 - 1, chosen_year)
    }
    df <- data.frame(
        sex=chosen_sex,
        year=chosen_year,
        N=1,
        after2005=as.numeric(chosen_year > 2005)
    ) %>% 
        mutate(after2005year=after2005*year)

    mu <- exp(predict(model_abx, newdata=df, type='link'))
    size <- exp(model_abx$family$getTheta())
    prob <- dnbinom(c(0:5), mu=mu, size=size)
    prob[6] <- 1 - sum(prob[1:5])
    return(data.frame(abx_exposure=c(0:5), prob_abx=prob))
}


OR_abx_calculator <- function(
    age, dose, params=c(1.711 + 0.115, -0.225, 0.053)
){
    if (dose == 0) {
        return(1)
    } else {
        return(exp(sum(params * c(1, pmin(age, 7), pmin(dose, 3)))))
    }
}


OR_fam_calculator <- function(
    age,
    fam_hist,
    params=c(log(OR_ASTHMA_AGE_3), (log(OR_ASTHMA_AGE_5) + log(OR_ASTHMA_AGE_3)) / 2)
){
    if (age < MIN_ASTHMA_AGE | fam_hist == 0 | age > 7) {
        return(1)
    } else {
        return(exp(params[1] + params[2] * (pmin(age, 5) - 3)))
    }
}


OR_risk_factor_calculator <- function(
    fam_hist,
    age,
    dose,
    params=list(
        c(log(OR_ASTHMA_AGE_3), (log(OR_ASTHMA_AGE_3) + log(OR_ASTHMA_AGE_5)) / 2 - log(OR_ASTHMA_AGE_3)),
        c(1.711 + 0.115, -0.225, 0.053)
    )
){
    if (age < MIN_ASTHMA_AGE) {
        return(1)
    } else {
        return(exp(
            log(OR_fam_calculator(age, fam_hist, params[[1]])) + 
            log(OR_abx_calculator(age, dose, params[[2]]))
        ))
    }
}


#' Compute the combined antibiotic exposure and family history odds ratio.
#' 
#' @param chosen_year The current year.
#' @param chosen_age The age of the person in years.
#' @param chosen_sex The sex of the person; 0 = female, 1 = male.
#' @param model_abx The fitted Negative Binomial model for the number of courses of antibiotics.
#' @param p_fam_distribution A dataframe with the probability of family history of asthma, given
#' that the person has asthma. Contains two columns: fam_history (0 or 1) and prob_fam.
#' @param df_fam_history_or A dataframe with the odds ratio of family history of asthma, given
#' the age of the person. Contains three columns: age (3, 4, or 5), fam_history (0 or 1),
#' and OR_fam: odds ratio.
#' @param df_abx_or A dataframe with the odds ratio of antibiotic exposure, given
#' the age of the person. Contains three columns: age (3, 4, or 5),
#' abx_exposure (0, 1, 2, 3, 4, or 5), and OR_abx: odds ratio.
#' @returns A dataframe with the following columns:
#' - fam_history: 0 or 1; 0 = no family history of asthma, 1 = family history of asthma
#' - abx_exposure: 0, 1, 2, 3, 4, 5(+); number of courses of antibiotics in the first year of life
#' - year: the current year
#' - sex: 0 or 1; 0 = female, 1 = male
#' - age: the age of the person in years
#' - prob: the probability of antibiotic exposure * probability of family history
#' - OR: the odds ratio of antibiotic exposure * odds ratio of family history
risk_factor_generator <- function(
    chosen_year, chosen_sex, chosen_age, model_abx, p_fam_distribution, df_fam_history_or, df_abx_or
){

    birth_year <- chosen_year - chosen_age
    df_abx_exposure <- p_antibiotic_exposure(max(birth_year, 2000), chosen_sex, model_abx)

    # combine abx_exposure = 3, 4, 5+ into 3+
    df_abx_exposure$prob_abx[4] <- sum(df_abx_exposure$prob_abx[4:6])
    df_abx_exposure <- df_abx_exposure %>% 
        filter(abx_exposure<=3)

    # select the given age if <= 5, otherwise select age == 5
    df_fam_history_or_age <- df_fam_history_or %>% 
        filter(age==min(chosen_age, 5)) %>% 
        select(-age)

    # select the given age if <= 8, otherwise select age == 8
    # filter out abx_exposure > 3
    df_abx_or_age <- df_abx_or %>% 
        filter(age == min(chosen_age, 8)) %>% 
        select(-age) %>%
        filter(abx_exposure<=3)        


    risk_set <- expand.grid(
        fam_history=c(0, 1),
        abx_exposure=c(0, 1, 2, 3, 4, 5)
    ) %>%
        mutate(
            year=chosen_year,
            sex=chosen_sex,
            age=chosen_age
        ) %>%
        filter(abx_exposure <= 3) %>% 
        left_join(p_fam_distribution, by=c("fam_history")) %>%
        left_join(df_abx_exposure, by=c("abx_exposure")) %>%
        left_join(df_fam_history_or_age, by=c("fam_history")) %>%
        left_join(df_abx_or_age, by=c("abx_exposure")) %>% 
        mutate(
            prob=prob_fam * prob_abx,
            OR=OR_abx * OR_fam
        ) %>%
        select(fam_history, abx_exposure, year, sex, age, prob, OR)
    return(risk_set)
}



# .5989652 -0.3574636

#' Compute the loss function given the effects of risk factors in the incidence equation, for
#' each year, age, and sex.
#' 
#' @param chosen_year The current year.
#' @param chosen_age The age of the person in years.
#' @param chosen_sex The sex of the person; 0 = female, 1 = male.
#' @param model_abx The fitted Negative Binomial model for the number of courses of antibiotics.
#' @param p_fam_distribution A dataframe with the probability of family history of asthma, given
#' that the person has asthma. Contains two columns: fam_history (0 or 1) and prob_fam.
#' @param df_fam_history_or A dataframe with the odds ratio of family history of asthma, given
#' the age of the person. Contains three columns: age (3, 4, or 5), fam_history (0 or 1),
#' and OR_fam: odds ratio.
#' @param df_abx_or A dataframe with the odds ratio of antibiotic exposure, given
#' the age of the person. Contains three columns: age (3, 4, or 5),
#' abx_exposure (0, 1, 2, 3, 4, or 5), and OR_abx: odds ratio.
#' @param df_incidence A dataframe with the incidence of asthma, with the following columns:
#' - year: the year
#' - age: the age in years
#' - sex: 0 or 1; 0 = female, 1 = male
#' - inc: the incidence of asthma
#' @param df_prevalence A dataframe with the prevalence of asthma, with the following columns:
#' - year: the year
#' - age: the age in years
#' - sex: 0 or 1; 0 = female, 1 = male
#' - prev: the prevalence of asthma
#' @param df_reassessment A dataframe with the reassessment of asthma, with the following columns:
#' - year: the year
#' - age: the age in years
#' - sex: 0 or 1; 0 = female, 1 = male
#' - ra: the reassessment of asthma
#' @param inc_beta_params A list of parameters for the incidence equation.
#' @param inc_function A function to compute the incidence correction.
#' @returns A list with the following elements:
#' - risk_set: a dataframe with the risk factors and their probabilities and odds ratios
#' - prev_sol: the solution for the prevalence calibration
#' - inc_sol: the solution for the incidence calibration
calibrator <- function(
    chosen_year,
    chosen_sex,
    chosen_age,
    model_abx,
    p_fam_distribution,
    df_fam_history_or,
    df_abx_or,
    df_incidence,
    df_prevalence,
    df_reassessment,
    inc_beta_params=c(0.3766256, -0.225),
    inc_function=inc_loss_function
){
  
    if(!is.list(inc_beta_params)){
        inc_beta_params <- list(
            c(log(OR_ASTHMA_AGE_3), inc_beta_params[1]),
            c(1.826, inc_beta_params[2], 0.053)
        )
    }

    risk_set <- risk_factor_generator(
        chosen_year, chosen_sex, chosen_age, model_abx, p_fam_distribution,
        df_fam_history_or, df_abx_or
    )

    target_prev <- df_prevalence %>% 
        filter(
            age==chosen_age & 
            year==chosen_year &
            sex==chosen_sex
        ) %>% 
        select(prev) %>% 
        unlist()
  
    if(chosen_age <= 7){

        prev_sol <- prev_calibrator(target_prev, risk_set$OR, risk_set$prob)
        p0 <- inverse_logit(logit(target_prev) - sum(risk_set$prob[-1] * prev_sol))
        risk_set$calibrated_prev <- inverse_logit(logit(p0) + log(risk_set$OR))
        risk_set$prev <- target_prev
        
        if(chosen_year== 2000){
            return(list(
                risk_set=risk_set,
                prev_sol=prev_sol,
                inc_sol=c()
            ))
        }
    
        if(chosen_age==3) {
            risk_set$calibrated_inc <- risk_set$calibrated_prev
            return(list(
                risk_set=risk_set,
                prev_sol=prev_sol,
                inc_sol=c()
            ))
        } else { # aged 4 or more
        
            risk_set$inc <- df_incidence %>% 
                filter(
                    age==chosen_age & 
                    year==chosen_year &
                    sex==chosen_sex
                ) %>% 
                select(inc) %>% 
                unlist()
      
            past_target_prev <- df_prevalence %>% 
                filter(
                    age==chosen_age - 1 & 
                    year==max(min_cal_year, chosen_year - 1) &
                    sex==chosen_sex
                ) %>% 
                select(prev) %>% 
                unlist()
        
            past_risk_set <- risk_factor_generator(
                max(min_cal_year, chosen_year - 1),
                chosen_sex,
                chosen_age - 1,
                model_abx,
                p_fam_distribution,
                df_fam_history_or,
                df_abx_or
            )
        
            target_RA <- df_reassessment %>% 
                filter(
                    age==chosen_age & 
                    year==chosen_year &
                    sex==chosen_sex
                ) %>% 
                select(ra) %>% 
                unlist()

            inc_risk_set <- risk_factor_generator(
                chosen_year, chosen_sex, chosen_age, model_abx, p_fam_distribution,
                df_fam_history_or, df_abx_or
            ) %>% 
                select(fam_history, abx_exposure, year, sex, age, prob)
        
            inc_risk_set$prob <- inc_risk_set$prob / sum(inc_risk_set$prob)
            
            inc_risk_set$OR <- inc_risk_set %>% 
                apply(., 1, FUN=function(x){
                    OR_risk_factor_calculator(
                        fam_hist=x[1],
                        age=x[5],
                        dose=x[2],
                        params=inc_beta_params
                    )
                })
        } 
    } else { # age > 7 => OR =1 for all abx
    
        risk_set <- risk_set %>% 
            group_by(fam_history, year, sex, age) %>% 
            summarise(
                prob=sum(prob),
                OR=mean(OR)
            ) %>% 
            ungroup()

        prev_sol <- prev_calibrator(target_prev, risk_set$OR, risk_set$prob)
        p0 <- inverse_logit(logit(target_prev) - sum(risk_set$prob[-1]*prev_sol))
        risk_set$calibrated_prev <- inverse_logit(logit(p0) + log(risk_set$OR))
        risk_set$prev <- target_prev
    
        if(chosen_year == 2000) {
            return(list(
                risk_set=risk_set,
                prev_sol=prev_sol,
                inc_sol=c()
            ))
        } else { 
 
            risk_set$inc <- df_incidence %>% 
                filter(
                    age==chosen_age & 
                    year==chosen_year &
                    sex==chosen_sex
                ) %>% 
                select(inc) %>% 
                unlist()
      
            past_target_prev <- df_prevalence %>% 
                filter(
                    age==chosen_age - 1 & 
                    year==max(min_cal_year, chosen_year - 1) &
                    sex==chosen_sex
                ) %>% 
                select(prev) %>% 
                unlist()

            past_risk_set <- risk_factor_generator(
                max(min_cal_year, chosen_year-1),
                chosen_sex,
                chosen_age - 1,
                model_abx,
                p_fam_distribution,
                df_fam_history_or,
                df_abx_or
            )

            if(chosen_age != 8){
                past_risk_set <- past_risk_set %>% 
                    group_by(fam_history) %>% 
                    summarise(
                        prob=sum(prob),
                        OR=mean(OR)
                    )
            } else {

                ttt_prev_sol <- prev_calibrator(past_target_prev, past_risk_set$OR, past_risk_set$prob)
                p0 <- inverse_logit(logit(past_target_prev) - sum(past_risk_set$prob[-1]*ttt_prev_sol))
                past_risk_set$calibrated_prev <- inverse_logit(logit(p0) + log(past_risk_set$OR))
                tmp_look <- past_risk_set %>% 
                    mutate(
                        yes_asthma=calibrated_prev * prob,
                        no_asthma=(1 - calibrated_prev) * prob
                    )
                past_tmp_OR <- (
                    sum(tmp_look$no_asthma[tmp_look$fam_history==0]) * 
                    sum(tmp_look$yes_asthma[tmp_look$fam_history==1]) /
                    (sum(tmp_look$yes_asthma[tmp_look$fam_history==0]) * 
                    sum(tmp_look$no_asthma[tmp_look$fam_history==1]))
                )
                past_risk_set <- past_risk_set %>% 
                    group_by(fam_history) %>% 
                    summarise(prob=sum(prob))
                past_risk_set$OR <- c(1, past_tmp_OR)
            }
            
            target_RA <- df_reassessment %>% 
                filter(
                    age==chosen_age & 
                    year==chosen_year &
                    sex==chosen_sex
                ) %>% 
                select(ra) %>% 
                unlist()
            
            inc_risk_set <- risk_factor_generator(
                chosen_year,chosen_sex,chosen_age, model_abx, p_fam_distribution, df_fam_history_or, df_abx_or
                ) %>% 
                filter(abx_exposure==0) %>% 
                select(fam_history,abx_exposure,year,sex,age,prob)
      
            inc_risk_set$OR <- inc_risk_set %>% 
                apply(., 1, FUN=function(x){
                    OR_risk_factor_calculator(
                        fam_hist=x[1],
                        age=x[5],
                        dose=x[2],
                        params=inc_beta_params
                    )
                })
        }
    }

    args <- list(
        target_inc=risk_set$inc[1],
        past_target_prev=past_target_prev,
        past_target_OR=past_risk_set$OR,
        target_OR=risk_set$OR,
        p_risk=past_risk_set$prob,
        ra=target_RA,
        misDx=0, # target misdiagnosis
        Dx=1, # target diagnosis
        risk_set=inc_risk_set,
        log_inc_OR=log(inc_risk_set$OR)[-1]
    )

    inc_sol <- do.call(inc_function, args)
  
    return(list(
        risk_set=risk_set,
        prev_sol=prev_sol,
        inc_sol=inc_sol
    ))
}


calculate_correction <- function(
    chosen_year,
    chosen_sex,
    chosen_age,
    model_abx,
    p_fam_distribution,
    df_fam_history_or,
    df_abx_or,
    df_incidence,
    df_prevalence,
    df_reassessment,
    inc_beta_params=optimized_inc_beta,
    inc_function=inc_correction_calculator
){

    df_results <- data.frame(
        year=chosen_year,
        sex=chosen_sex,
        age=chosen_age,
        obj_value=NA,
        prev_correction=NA,
        inc_correction=NA
    )
    results <- calibrator(
        chosen_year,
        chosen_sex,
        chosen_age,
        model_abx,
        p_fam_distribution,
        df_fam_history_or,
        df_abx_or,
        df_incidence,
        df_prevalence,
        df_reassessment,
        inc_beta_params=inc_beta_params,
        inc_function=inc_function
    )
    risk_set <- results$risk_set
    prev_sol <- results$prev_sol
    inc_sol <- results$inc_sol
        df_results$prev_correction <- -sum(risk_set$prob[-1]*prev_sol)
    if (chosen_year > 2000) {
    if(chosen_age==3) {
        df_results$inc_correction <- -sum(risk_set$prob[-1]*prev_sol)
    } else { # aged 4 or more
        df_results$obj_value <- inc_sol[1]
        df_results$inc_correction <- inc_sol[2]
    }
    }
    return(df_results)
}


generate_correction <- function(
    df,
    inc_beta_params,
    model_abx,
    p_fam_distribution,
    df_fam_history_or,
    df_abx_or,
    df_incidence,
    df_prevalence,
    df_reassessment
){
    apply(df, 1, FUN=function(x) {
        calculate_correction(
            chosen_year=x[1],
            chosen_sex=x[2],
            chosen_age=x[3],
            df_incidence=df_incidence,
            df_prevalence=df_prevalence,
            df_reassessment=df_reassessment,
            p_fam_distribution=p_fam_distribution,
            df_fam_history_or=df_fam_history_or,
            df_abx_or=df_abx_or,
            model_abx=model_abx,
            inc_beta_params=inc_beta_params
        )
    }) %>% 
    do.call(rbind,.)
}


inc_beta_solver <- function(
    model_abx,
    df_fam_history_or,
    df_abx_or,
    df_incidence,
    df_prevalence,
    df_reassessment,
    baseline_year=2001,
    stabilization_year=2025,
    max_age=63,
    inc_beta_params=INC_BETA_PARAMS
){
    cal_years <- baseline_year:(stabilization_year+1)
    ages <- 4:max_age
    sexes <- 0:1
    covar <- expand.grid(year=cal_years,sex=sexes,age=ages) %>% 
        as.data.frame()

    obj <- function(inc_beta_params){
        apply(covar, 1, FUN=function(x) {
            calibrator(
                x[1], x[2], x[3], model_abx, p_fam_distribution, df_fam_history_or,
                df_abx_or, df_incidence, df_prevalence, df_reassessment, inc_beta_params
            )$inc_sol
        }) %>% mean()
    }
  
    res_optim <- optim(unlist(inc_beta_params),fn=obj,method='BFGS')
    res_nlm <- nlm(obj,unlist(inc_beta_params),steptol=1e-6,gradtol=1e-6,print.level=2)
    write_rds(res_optim, here("R/res_optim.rds"))
}


# asthma prev and inc -----------------------------------------------------

asthma_inc_model <- read_rds(here("R/asthma_incidence_model.rds"))
asthma_prev_model <- read_rds(here("R/asthma_prevalence_model.rds"))


df_asthma <- expand.grid(
    age=3:110,
    sex=c(0, 1),
    year=min_cal_year:max_cal_year
) %>% 
    as.data.frame()

df_asthma <- df_asthma %>% 
    mutate(inc=asthma_predictor(age, sex, year, "inc")) %>% 
    mutate(prev=asthma_predictor(age, sex, year, "prev")) %>% 
    mutate(inc=ifelse(age==3, prev, inc))

df_incidence <- df_asthma %>% 
    select(year, age, sex, inc) %>% 
    pivot_wider(names_from=sex, values_from=inc) %>% 
    as.data.frame()
colnames(df_incidence)[c(3, 4)] <- c("F", "M")
df_incidence$province <- chosen_province

df_incidence <- df_incidence %>% 
    select(-province) %>%
    pivot_longer(3:4, values_to="inc", names_to='sex') %>% 
    mutate(sex=as.numeric(sex=="M"))

df_prevalence <- df_asthma %>% 
    select(year, age, sex, prev) %>% 
    pivot_wider(names_from=sex, values_from=prev) %>% 
    as.data.frame()
colnames(df_prevalence)[c(3, 4)] <- c("F", "M")
df_prevalence$province <- chosen_province

df_prevalence <- df_prevalence %>% 
    select(-province)%>% 
    pivot_longer(3:4, values_to="prev", names_to='sex')%>% 
    mutate(sex=as.numeric(sex=="M"))

df_reassessment <- read_csv(here("src/processed_data/master_asthma_reassessment.csv")) %>% 
    filter(province==chosen_province)

df_reassessment <- df_reassessment %>% 
    select(-province)%>% 
    pivot_longer(3:4, values_to="ra", names_to='sex')%>% 
    mutate(sex=as.numeric(sex=="M"))


# risk factors ------------------------------------------------------------

p_fam_distribution <- data.frame(
    fam_history=c(0, 1),
    prob_fam=c(1 - PROB_FAM_HIST, PROB_FAM_HIST)
)

# Abx exposure: 0 1 2 3 4 5+
# differs by year
model_abx <- read_rds(here("R/BC_count_model.rds"))

# prev eqn OR
df_fam_history_or <- list(
    c(1, OR_ASTHMA_AGE_3),
    c(1, exp((log(OR_ASTHMA_AGE_3) + log(OR_ASTHMA_AGE_5)) / 2)),
    c(1, OR_ASTHMA_AGE_5)
)
df_fam_history_or <- data.frame(
    age=c(3, 4, 5), do.call(rbind, df_fam_history_or)
)
colnames(df_fam_history_or)[-1] <- c(0,1)
df_fam_history_or <- pivot_longer(
    df_fam_history_or, cols=-1, names_to="fam_history", values_to="OR_fam"
) %>% 
    mutate(fam_history=as.numeric(fam_history))

# fam_history + \beta_age * (age-3) + dose()
# free parameters are : age

df_abx_or <- read_csv(here("R/dose_response_log_aOR.csv"))
colnames(df_abx_or) <- c("age", paste0("OR", c(1:5)))
df_abx_or$OR0 <- 0
df_abx_or <- df_abx_or %>% 
    select(age, OR0, OR1:OR5) %>% 
    mutate(across(contains("OR"), exp)) %>% 
    filter(age >= 3) %>% 
    rbind(data.frame(age=8, OR0=1, OR1=1, OR2=1, OR3=1, OR4=1, OR5=1))

df_abx_or <- pivot_longer(
    df_abx_or, cols=-1, names_to="abx_exposure", values_to="OR_abx"
) %>% 
    mutate(abx_exposure=as.numeric(str_remove(abx_exposure, "OR")))

# incorporate the estimates of the risk factors and correction terms ----------
baseline_year=2001
stabilization_year=2025
max_age=63

inc_beta_solver(
    df_incidence,
    df_prevalence,
    df_reassessment,
    p_fam_distribution,
    df_fam_history,
    df_abx,
    model_abx,
    baseline_year=baseline_year,
    stabilization_year=stabilization_year,
    max_age=max_age
)

res_optim <- read_rds(here("R/res_optim.rds"))
optimized_inc_beta <- res_optim$par

cal_years <- (baseline_year-1):(stabilization_year+1)
ages <- 3:max_age
sexes <- 0:1
calibration_results <- expand.grid(year=cal_years,sex=sexes,age=ages) %>%
  as.data.frame()


df_correct <- generate_correction(calibration_results, model_abx, optimized_inc_beta)

df_correct_prev <- df_correct %>% 
    select(1:3, 5) %>% 
    rename(correction=prev_correction) %>% 
    mutate(type='prev')

df_correct_inc <- df_correct %>% 
    select(1:3, 6) %>% 
    rename(correction=inc_correction) %>% 
    mutate(type='inc')

master_correct <- rbind(df_correct_prev, df_correct_inc)
master_correct <- master_correct %>% 
  mutate(correction=ifelse(is.na(correction), 0, correction))

write_csv(master_correct, here("src/processed_data/master_asthma_occurrence_correction.csv"))

# examine_results <- df_correct %>% 
#   filter(age !=3)
# 
# ggplot(examine_results %>%
#          filter(age <= 10) %>% 
#          mutate(year=as.factor(year)),aes(y=obj_value,x=age,col=year)) +
#   geom_point() + 
#   facet_grid(.~sex)
# # write_csv(final_result,"master_calibrated_asthma_prev_inc_M3.csv")
# 
# #### post-processing
# final_result <- read_csv("master_calibrated_asthma_prev_inc_M3.csv")
# 
# tmp <- head(final_result %>% filter(year==2004),n=8)
# tmp$calibrated_prev[seq(2,8,by=2)]/tmp$calibrated_prev[seq(1,7,by=2)]
# 
# logit <- function(p){
#   log(p/(1-p))
# }
# 
# inv_logit <- function(x){
#   exp(x)/(1+exp(x))
# }
# 
# final_result %>% 
#   group_by(fam_history,abx_exposure,sex,age) %>% 
#   summarise(mean(OR),median(OR),sd(OR)) %>% 
#   filter(`sd(OR)`!=0)  # variations are negligible 
# 
# tmp <- final_result %>% 
#   filter(age==5 & sex==0 & fam_history==0 & abx_exposure==0)
# 