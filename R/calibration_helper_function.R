library(tidyverse)
library(epitools)
library(mgcv)

logit <- function(p){
  log(p/(1-p))
}

inverse_logit <-function(x){
  exp(x)/(1+exp(x))
}

OR_generator <- function(risk_set,params){
  risk_set$OR <- apply(risk_set,1,FUN = function(x){
    exp(sum(x*params))
  })
  return(risk_set)
}


#' @title compute_asthma_prev_risk_factors
#' @description This function calculates the asthma prevalence based on the risk factors and the
#' parameters provided by x.
#' @param x A vector of parameters to be optimized.
#' @param multiple_risk_factors A boolean indicating if there are multiple risk factors.
#' @param target_OR A vector of odds ratios between the risk factors and asthma.
#' @param risk_factor_prev A vector of the prevalence of the risk factor levels.
#' @param beta0 The intercept of the logistic regression model.
#' @return The calibrated asthma prevalence.
#' @details This function is used internally by the prev_calibrator function.
compute_asthma_prev_risk_factors <- function(
    x, multiple_risk_factors, target_OR, risk_factor_prev, beta0
) {
    # binary
        if(!multiple_risk_factors) {
            if(length(target_OR) == 1) {
            p0 <- inverse_logit(beta0 - risk_factor_prev * x)
            asthma_prev_x <- inverse_logit(logit(p0) + log(target_OR))
            return(
                sum(asthma_prev_x * risk_factor_prev) + 
                (1 - risk_factor_prev) * p0
            )
            } else {
            # asthma_prev_x: asthma prevalence at risk factor level x
            asthma_prev_x <- inverse_logit(beta0 + log(target_OR) - sum(risk_factor_prev[-1] * x))
            return(sum(asthma_prev_x * risk_factor_prev))
            }
        } else {
      # number of risk factors
        n_risk_factors <- length(risk_factor_prev)
      # number of levels of risk factors
        p_length <- lapply(risk_factor_prev, length) %>% unlist()
      p_length_optim <- p_length - 1 
      
      # break up x
        indices <- c()
        xs <- c()
      index_start <- 1
        for(i in 1:n_risk_factors){
        index_end <- index_start + p_length_optim[[i]] - 1 
        xs[[i]] <- x[index_start:index_end]
                index_start <- index_end + 1
      } 
      
            penalty <- mapply(
            function(tmp_risk_factor_prev, tmp_x){
                sum(tmp_risk_factor_prev[-1] * tmp_x)
                },
            risk_factor_prev,
                xs,
                SIMPLIFY=FALSE
            ) %>% unlist()
      
        asthma_prev_x_unlisted <- inverse_logit(beta0 - sum(penalty) + log(unlist(target_OR)))
        risk_factor_prev_unlisted <- unlist(risk_factor_prev)
        return(sum(asthma_prev_x_unlisted * risk_factor_prev_unlisted))
    }
}

#' @title compute_asthma_prevalence_difference
#' @description This function calculates the objective function for the optimization process.
#' @param x A vector of parameters to be optimized.
#' @param multiple_risk_factors A boolean indicating if there are multiple risk factors.
#' @param target_OR A vector of odds ratios between the risk factors and asthma.
#' @param risk_factor_prev A vector of the prevalence of the risk factor levels.
#' @param beta0 The intercept of the logistic regression model.
#' @param asthma_prev_target The target prevalence of asthma.
#' @return The absolute difference between the calculated and target prevalence.
#' @details This function is used internally by the prev_calibrator function.
compute_asthma_prevalence_difference <- function(
    x, multiple_risk_factors, target_OR, risk_factor_prev, beta0, asthma_prev_target
) {
    
    asthma_prev_calibrated <- compute_asthma_prev_risk_factors(
        x, multiple_risk_factors, target_OR, risk_factor_prev, beta0
    )
    return(abs(asthma_prev_calibrated - asthma_prev_target))
    }


#' @title prev_calibrator
#' @description This function calibrates the prevalence of asthma in a population
#'   based on the target prevalence and odds ratios of risk factors.
#' @param asthma_prev_target The target prevalence of asthma.
#' @param target_OR A vector of odds ratios for the risk factors.
#' @param risk_factor_prev A vector of the prevalence of the risk factors.
#' @param beta0 The intercept of the logistic regression model.
#' @param multiple_risk_factors A boolean indicating if there are multiple risk factors.
#' @param chosen_trace A boolean indicating if the trace should be printed.
#' @return A vector of the calibrated asthma prevalence for each risk factor level.
prev_calibrator <- function(
    asthma_prev_target,
    target_OR,
    risk_factor_prev,
    beta0=NULL,
    multiple_risk_factors=FALSE,
    chosen_trace=FALSE
) {

    if(is.null(beta0)){
        beta0 <- logit(asthma_prev_target)
    }

  if(!multiple_risk_factors){
        if (length(target_OR)==1){
            n_params <- 1
        } else {
            n_params <- length(target_OR) - 1
        }
    } else {
        n_params <- sum(unlist(target_OR) != 1)
  }

    return(optim(
        par=rep(0, n_params),
        fn=compute_asthma_prevalence_difference,
        multiple_risk_factors=multiple_risk_factors,
        target_OR=target_OR,
        risk_factor_prev=risk_factor_prev,
        beta0=beta0,
        asthma_prev_target=asthma_prev_target,
        control=list(abstol=1e-15, maxit=10000, trace=chosen_trace),
        method="BFGS",
        hessian=TRUE)$par
    ) 
}


# return 1) correction term for prev
#        2) correction term for inc
#        3) OR for inc


#' @title compute_contingency_table
#' @description This function generates a table of the proportions of the population
#' at different levels of family history and antibiotic exposure.
#' @param risk_factor_prev A vector of the prevalence of the risk factor levels.
#' @param target_OR A vector of odds ratios for the risk factors.
#' @param asthma_prev_calibrated A vector of the calibrated asthma prevalence.
#' @return A list of vectors representing the proportions of the population for different risk
#' factor levels / combinations. For example, if we have the risk factors of family history {0, 1}
#' and antibiotic exposure {0, 1, 2, 3}, then we have 2 * 4 = 8 combinations. Each combination is
#' called a "risk factor level" and is indexed by i (this corresponds to the index in the risk_set
#' table). The first combination, i = 1 is a special case; this is where there are no risk factors.
#' We use this combination, referred to as the "ref" level, in the calculation of all the tables.
#' Each list entry contains a vector of length 4, with the following entries:
#' - di: proportion of population labelled as no asthma with no risk factors
#' - ci: proportion of population labelled as asthma with no risk factors
#' - bi: proportion of population labelled as no asthma with risk factors
#' - ai: proportion of population labelled as asthma with risk factors
#' @details This function is used internally by the inc_correction_calculator function.
#' @note The function uses the metafor package to convert odds ratios into proportions.
compute_contingency_table <- function(
    risk_factor_prev, target_OR, asthma_prev_calibrated
) {
    prev_table <- c()

    asthma_prev_ref <- asthma_prev_calibrated[1]
    risk_factor_prev_ref <- risk_factor_prev[1]

    for(i in 2:(length(target_OR))){

        # prevalence of risk factor combination i
        risk_factor_prev_i <- risk_factor_prev[i] / (risk_factor_prev[i] + risk_factor_prev_ref)
        # calibrated asthma prevalence
        asthma_prev <- asthma_prev_calibrated[i]

        # return: a b c d
        # Solve the following:
        # tmp_target_prev  = (b+d)/(a+b+c+d) 
        # tmp_p[1] = b/(a+b) 
        # tmp_p[2] = d/(c+d) 
        # tmp_target_OR  = (a*d)/(b*c) 
    
        # metafor pkg
        # Bonett, D. G. (2007).
        # Transforming odds ratios into correlations for meta-analytic research. 
        # American Psychologist, 62(3), 254–255. ⁠https://doi.org/10.1037/0003-066x.62.3.254⁠

        #                | asthma | no asthma |
        # --------------------------------------------
        # risk factor    |   ai   |     bi    |  n1i
        # --------------------------------------------
        # no risk factor |   ci   |     di    |
        # --------------------------------------------
        #                |  n2i   |           |  ni


        sample_size <- 1e10
        prev_table[[i]] <- rev(
            metafor::conv.2x2(
                ori=target_OR[i],
                ni=sample_size,
                n1i=risk_factor_prev_i * sample_size, # prev of exposure
                n2i=sum(c(1 - risk_factor_prev_i, risk_factor_prev_i) * c(asthma_prev_ref, asthma_prev)) * sample_size # prev of asthma
            ) / sample_size
        )
    }
    return(prev_table)
}


#' @title compute_odds_ratio_difference
#' @description This function calculates the difference in odds ratios between the target and the
#'   calibrated asthma prevalence.
#' @param asthma_inc_calibrated A vector of the calibrated asthma incidence.
#' @param target_OR A vector of odds ratios for the risk factors.
#' @param contingency_table_past A list of vectors representing the proportions of the population
#' for different risk factor levels / combinations for the previous year.
#' For example, if we have the risk factors of family history {0, 1} and 
#' antibiotic exposure {0, 1, 2, 3}, then we have 2 * 4 = 8 combinations. Each combination is
#' called a "risk factor level" and is indexed by i (this corresponds to the index in the risk_set
#' table). The first combination, i = 1 is a special case; this is where there are no risk factors.
#' We use this combination, referred to as the "ref" level, in the calculation of all the tables.
#' Each list entry contains a vector of length 4, with the following entries:
#' - di: proportion of population labelled as no asthma with no risk factors
#' - ci: proportion of population labelled as asthma with no risk factors
#' - bi: proportion of population labelled as no asthma with risk factors
#' - ai: proportion of population labelled as asthma with risk factors
#' @param ra_target A value between 0 and 1 indicating the target reassessment.
#' @param misDx A numeric value representing the misdiagnosis rate.
#' @param Dx A numeric value representing the diagnosis rate.
#' @return The mean difference in log odds ratios.
#' @details This function is used internally by the inc_correction_calculator function.
compute_odds_ratio_difference <- function(
    asthma_inc_calibrated,
    target_OR,
    contingency_table_past,
    ra_target=1,
                           misDx=0,
    Dx=1
    ) {

    asthma_inc_calibrated_ref <- asthma_inc_calibrated[1]
    total_diff_log_OR <- 0
    
    for(i in 2:(length(target_OR))){
        asthma_inc <- asthma_inc_calibrated[i]
      
        # contingency table of the population with asthma from a previous year
        ref_a0 <- contingency_table_past[[i]][1] # proportion of population labelled as no asthma with no risk factors
        ref_b0 <- contingency_table_past[[i]][2] # proportion of population labelled as asthma with no risk factors
        ref_c0 <- contingency_table_past[[i]][3] # proportion of population labelled as no asthma with risk factors level i
        ref_d0 <- contingency_table_past[[i]][4] # proportion of population labelled as asthma with risk factors level i
      
      # contingency table of the population with asthma from a previous year
      # if ra=1, no reversibility
        a0 <- ref_b0 * (1 - ra_target) # proportion of population who lose asthma diagnosis with no risk factors
        c0 <- ref_d0 * (1 - ra_target) # proportion of population who lose asthma diagnosis at risk factors level i
        b0 <- ref_b0 * ra_target # proportion of population who keep asthma diagnosis with no risk factors
        d0 <- ref_d0 * ra_target # proportion of population who keep asthma diagnosis at risk factors level i
      
      # contingency table of the exposure level 
        # no risk factors & no asthma: 
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * correct Dx = no asthma + 
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * misDx = no asthma
        a1 <- ref_a0 * ((1 - asthma_inc_calibrated_ref) * (1 - misDx) + asthma_inc_calibrated_ref * (1 - Dx))
        # no risk factors & yes asthma: get asthma and correctly Dx + did not get asthma but misDx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * misdiagnosis = has asthma + 
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * correct Dx = has asthma
        b1 <- ref_a0 * ((1 - asthma_inc_calibrated_ref) *  misDx + asthma_inc_calibrated_ref * Dx)
        # yes risk factors & no asthma: did not get asthma and correctly Dx + got asthma but incorrectly Dx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * correct Dx = no asthma +
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * misDx = no asthma
        c1 <- ref_c0 * ((1 - asthma_inc) * (1 - misDx) + asthma_inc * (1 - Dx))
        # yes risk factors & yes asthma: got asthma and correctly Dx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * misDx = has asthma +
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * correct Dx = has asthma
        d1 <- ref_c0 * ((1 - asthma_inc) * misDx + asthma_inc * Dx)
      
      # two targets
            # objective: asthma prev OR
        # (no risk factors) proportion of population who either: 
        # (a0) lose asthma diagnosis from t0 - t1 or (a1) do not get new asthma diagnosis at t1
      a <- a0 + a1
        # (no risk factors) proportion of population who either:
        # (b0) keep asthma diagnosis from t0 - t1 or (b1) get new asthma diagnosis at t1
      b <- b0 + b1
        # (risk factors i) proportion of population who either:
        # (c0) lose asthma diagnosis from t0 - t1 or (c1) do not get new asthma diagnosis at t1
      c <- c0 + c1
        # (risk factors i) proportion of population who either:
        # (d0) keep asthma diagnosis from t0 - t1 or (d1) get new asthma diagnosis at t1
      d <- d0 + d1

        # odds ratio = (a*d)/(b*c)
        diff_log_OR <- abs(log(target_OR[i]) - (log(d) + log(a) - log(b) - log(c)))
        total_diff_log_OR <- total_diff_log_OR + diff_log_OR
    }
    
    return(total_diff_log_OR %>% unlist() %>% mean())
}

inc_correction_calculator <- function(
    asthma_inc_target,
    asthma_prev_target_past,
    past_target_OR,
    target_OR,
    risk_factor_prev_past,
    risk_set,
    ra_target=1.0,
    misDx=0,
    Dx=1,
){
  
    beta0 <- logit(asthma_inc_target)
    log_inc_OR <- log(risk_set$OR)
    
    # asthma prevalance ~ risk factor parameters for the previous year
    asthma_prev_risk_factor_params_past <- prev_calibrator(
        asthma_prev_target=asthma_prev_target_past,
        target_OR=past_target_OR,
        risk_factor_prev=risk_factor_prev_past
    )

    # calibrated asthma prevalence for the previous year
    asthma_prev_calibrated_past <- inverse_logit(
        logit(asthma_prev_target_past) +
        log(past_target_OR) - 
        sum(risk_factor_prev_past[-1] * asthma_prev_risk_factor_params_past) 
    )
    # distribution of the risk factors for the population without asthma
    risk_factor_prev_past_no_asthma <- (1 - asthma_prev_calibrated_past) * risk_factor_prev_past
    # normalize
    risk_factor_prev_past_no_asthma <- risk_factor_prev_past_no_asthma / sum(risk_factor_prev_past_no_asthma)
    
    # asthma prevalance ~ risk factor parameters for incidence
    asthma_prev_risk_factor_params <- prev_calibrator(
        asthma_prev_target=asthma_inc_target,
        target_OR=exp(log_inc_OR),
        risk_factor_prev=risk_factor_prev_past_no_asthma
    )

    asthma_inc_correction <- sum(
        asthma_prev_risk_factor_params * risk_factor_prev_past_no_asthma[-1]
    )

    # calibrated asthma incidence
    asthma_inc_calibrated <- inverse_logit(beta0 + log_inc_OR - asthma_inc_correction)
    
    # for each OR, we need to obtain the contingency table
    contingency_table <- compute_contingency_table(
        risk_factor_prev=risk_factor_prev_past,
        target_OR=past_target_OR,
        asthma_prev_calibrated=asthma_prev_calibrated_past
    )
    
    mean_diff_log_OR <- compute_odds_ratio_difference(
        asthma_inc_calibrated, target_OR, contingency_table, ra_target, misDx, Dx
    )

    return(list(
        mean_diff_log_OR=mean_diff_log_OR,
        asthma_inc_correction=-asthma_inc_correction
    ))
}