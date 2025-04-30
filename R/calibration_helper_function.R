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
#' @return A list of vectors representing the proportions of the population at different levels.
#' The list index corresponds to the index in the risk_set table. Each list entry contains a vector
#' of length 4, with the following entries:
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


  
obj_function <- function(
    log_inc_OR,
    no_asthma_risk_factor_prev_dist,
    target_OR,
    asthma_inc_target,
    beta0,
    prev_table,
                           ra=1,
                           misDx=0,
    Dx=1
    ) {

    # calibrate the current incidence to the target incidence
    asthma_prev_risk_factor_params <- prev_calibrator(
            asthma_prev_target=asthma_inc_target,
        target_OR=exp(log_inc_OR),
            risk_factor_prev=no_asthma_risk_factor_prev_dist
        )

    asthma_inc_correction <- sum(
        asthma_prev_risk_factor_params * no_asthma_risk_factor_prev_dist[-1]
    )

    asthma_inc_calibrated <- inverse_logit(beta0 + log_inc_OR - asthma_inc_correction)
    asthma_inc_calibrated_ref <- asthma_inc_calibrated[1]
    
    total_diff_log_OR <- 0
    
    for(i in 2:(length(target_OR))){
        asthma_inc <- asthma_inc_calibrated[i]
      
        ref_a0 <- prev_table[[i]][1] # proportion of population labelled as no asthma at reference level
        ref_b0 <- prev_table[[i]][2] # proportion of population labelled as asthma at reference level
        ref_c0 <- prev_table[[i]][3] # proportion of population labelled as no asthma at level x
        ref_d0 <- prev_table[[i]][4] # proportion of population labelled as asthma at level x
      
      # contingency table of the population with asthma from a previous year
      # if ra=1, no reversibility
        a0 <- ref_b0 * (1 - ra) # proportion of population who lose asthma diagnosis at reference level
        c0 <- ref_d0 * (1 - ra) # proportion of population who lose asthma diagnosis at level x
        b0 <- ref_b0 * ra # proportion of population who keep asthma diagnosis at reference level
        d0 <- ref_d0 * ra # proportion of population who keep asthma diagnosis at level x
      
      # contingency table of the exposure level 
        # no antibiotic exposure & no asthma at reference level: 
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * correct Dx = no asthma + 
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * misDx = no asthma
        a1 <- ref_a0 * ((1 - asthma_inc_calibrated_ref) * (1 - misDx) + asthma_inc_calibrated_ref * (1 - Dx))
        # no antibiotic exposure & yes asthma: get asthma and correctly Dx + did not get asthma but misDx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * misdiagnosis = has asthma + 
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * correct Dx = has asthma
        b1 <- ref_a0 * ((1 - asthma_inc_calibrated_ref) *  misDx + asthma_inc_calibrated_ref * Dx)
        # yes antibiotic exposure & no asthma: did not get asthma and correctly Dx + got asthma but incorrectly Dx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * correct Dx = no asthma +
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * misDx = no asthma
        c1 <- ref_c0 * ((1 - asthma_inc) * (1 - misDx) + asthma_inc * (1 - Dx))
      # yes exposure & yes asthma: got asthma and correctly Dx
        # t0 = no asthma diagnosis, t1 = no asthma diagnosis * misDx = has asthma +
        # t0 = no asthma diagnosis, t1 = asthma diagnosis * correct Dx = has asthma
        d1 <- ref_c0 * ((1 - asthma_inc) * misDx + asthma_inc * Dx)
      
      # two targets
            # objective: asthma prev OR
        # (ref) proportion of population who either lose asthma diagnosis from t0 - t1 or do not get new asthma diagnosis at t1
      a <- a0 + a1
        # (ref) proportion of population who either keep asthma diagnosis from t0 - t1 or get new asthma diagnosis at t1
      b <- b0 + b1
        # (level x) proportion of population who either lose asthma diagnosis from t0 - t1 or do not get new asthma diagnosis at t1
      c <- c0 + c1
        # (level x) proportion of population who either keep asthma diagnosis from t0 - t1 or get new asthma diagnosis at t1
      d <- d0 + d1

        # odds ratio = (a*d)/(b*c)
        print(paste0("a: ", a))
        print(paste0("b: ", b))
        print(paste0("c: ", c))
        print(paste0("d: ", d))
        print(target_OR[i])
        diff_log_OR <- abs(log(target_OR[i]) - (log(d) + log(a) - log(b) - log(c)))
        total_diff_log_OR <- total_diff_log_OR + diff_log_OR
        print(total_diff_log_OR)
    }
    
    return(
        list(
            mean_diff_log_OR=total_diff_log_OR %>% unlist() %>% mean(),
            asthma_inc_correction=-asthma_inc_correction
        )
    )
}

inc_correction_calculator <- function(
    asthma_inc_target,
    asthma_prev_target_past,
    past_target_OR,
    target_OR,
    risk_factor_prev,
    ra=1,
    misDx=0,
    Dx=1,
    risk_set
){
  
    beta0 <- logit(asthma_inc_target)
    log_inc_OR <- log(risk_set$OR)
    prop_asthma <- asthma_prev_target_past # proportion with asthma
    prop_no_asthma <- (1 - asthma_prev_target_past) # proportion without asthma
    
    # reconstruct contingency table for each OR
    asthma_prev_risk_factor_params <- prev_calibrator(
        asthma_prev_target=asthma_prev_target_past,
        target_OR=past_target_OR,
        risk_factor_prev=risk_factor_prev
    )

    asthma_prev_calibrated <- inverse_logit(
        logit(asthma_prev_target_past) +
        log(past_target_OR) - 
        sum(risk_factor_prev[-1] * asthma_prev_risk_factor_params) 
    )
    # distribution of the risk factors for the population without asthma
    no_asthma_risk_factor_prev_dist <- (1 - asthma_prev_calibrated) * risk_factor_prev
    # normalize
    no_asthma_risk_factor_prev_dist <- no_asthma_risk_factor_prev_dist / sum(no_asthma_risk_factor_prev_dist)
    
    # for each OR, we need to obtain the contingency table
    contingency_table <- compute_contingency_table(
        risk_factor_prev=risk_factor_prev,
        target_OR=past_target_OR,
        asthma_prev_calibrated=asthma_prev_calibrated
    )
    
    fnc_value <- obj_function(
        log_inc_OR, no_asthma_risk_factor_prev_dist, target_OR, asthma_inc_target, beta0,
        contingency_table, ra, misDx, Dx
    )
    return(fnc_value)
}