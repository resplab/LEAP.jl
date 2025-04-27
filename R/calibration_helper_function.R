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

inc_loss_function <- function(
    target_inc,
                           past_target_prev,
                           past_target_OR,
                           target_OR,
                           p_risk,
                           ra=1,
                           misDx=0,
                           Dx=1,
                           risk_set,
    log_inc_OR
){
  
  beta0 <- logit(target_inc)
    prop_asthma <- past_target_prev # proportion with asthma
    prop_no_asthma <- (1 - past_target_prev) # proportion without asthma
  
  # reconstruct contingency table for each OR
  ref_p_risk <- p_risk[1]
    sol <- prev_calibrator(
        asthma_prev_target=past_target_prev,
        target_OR=past_target_OR,
        risk_factor_prev=p_risk
    )
    p0 <- inverse_logit(logit(past_target_prev) - sum(p_risk[-1] * sol))
  calibrated_p <- inverse_logit(logit(p0) + log(past_target_OR))
  # distribution of the risk factors for the population without asthma
    no_asthma_p_risk_dist <- (1 - calibrated_p) * p_risk
  # normalize
    no_asthma_p_risk_dist <- no_asthma_p_risk_dist / sum(no_asthma_p_risk_dist)
  
  # for each OR, we need to obtain the contingency table
  
  prev_table <- c()
  
    for(i in 1:(length(past_target_OR) - 1)){
        tmp_p_risk <- p_risk[c(1, i + 1)]
        tmp_p_risk <- tmp_p_risk / sum(tmp_p_risk)
        tmp_p <- calibrated_p[c(1, i + 1)]

    # return: a b c d
    # a: no exp, no asthma
    # b: no exp, yes asthma
    # c: yes exp, no asthma
    # d: yes exp, yes asthma
    # Solve the following:
    # tmp_target_prev  = (b+d)/(a+b+c+d) 
    # tmp_p[1] = b/(a+b) 
    # tmp_p[2] = d/(c+d) 
    # tmp_target_OR  = (a*d)/(b*c) 
    
        # metafor pkg
    # Bonett, D. G. (2007).
    # Transforming odds ratios into correlations for meta-analytic research. 
    # American Psychologist, 62(3), 254–255. ⁠https://doi.org/10.1037/0003-066x.62.3.254⁠
    
    nn <- 1e10
        prev_table[[i]] <- rev(
            metafor::conv.2x2(
                ori=past_target_OR[i + 1],
                ni=nn,
                n1i=((1 - tmp_p[2]) * tmp_p_risk[2] + tmp_p_risk[2] * tmp_p[2]) * nn, # prev of exposure
                n2i=sum(tmp_p_risk * tmp_p) * nn # prev of asthma
            ) / nn
        )
    }

    target_prev <- (
        past_target_prev * ra + 
        target_inc * (1 - past_target_prev) * Dx +
        (1-target_inc)*(1-past_target_prev)*misDx
    )
    tmp_sol <- prev_calibrator(target_prev, target_OR, p_risk)
    tmp_p0 <- inverse_logit(logit(target_prev) - sum(p_risk[-1] * tmp_sol))
  tmp_calibrated_p <- inverse_logit(logit(tmp_p0) + log(target_OR))
  
    future_prev_table <- c()

    for(i in 1:(length(target_OR) - 1)){
        tmp_p_risk <- p_risk[c(1, i + 1)]
        tmp_p_risk <- tmp_p_risk / sum(tmp_p_risk)
        tmp_p <- tmp_calibrated_p[c(1, i + 1)]
    nn <- 1e10
        prev_table[[i]] <- rev(
            metafor::conv.2x2(
                ori=target_OR[i + 1],
                ni=nn,
                n1i=((1 - tmp_p[2]) * tmp_p_risk[2] + tmp_p_risk[2] * tmp_p[2]) * nn, # prev of exposure
                n2i=sum(tmp_p_risk * tmp_p) * nn # prev of asthma
            ) / nn
        )
  }
  
    obj_function <- function(
        y, no_asthma_p_risk_dist, target_OR, target_inc, beta0, prev_table
    ) {
        x <- y
        q <- length(no_asthma_p_risk_dist) - 1
    target_OR_no_ref <- target_OR[-1]
        # calibrate the current inc to the target inc
        tmp_sol <- prev_calibrator(
            asthma_prev_target=target_inc,
            target_OR=exp(c(0, x)),
            risk_factor_prev=no_asthma_p_risk_dist
        )
        logit_p0 <- beta0 - sum(tmp_sol * no_asthma_p_risk_dist[-1])
        calibrated_inc <- inverse_logit(logit_p0 + c(0, x))
    
    ref_cal_inc <- calibrated_inc[1]
    calibrated_inc_no_ref <- calibrated_inc[-1]
    
    result <- 0
    
        for(i in 1:(length(target_OR) - 1)){
      cal_inc <- calibrated_inc_no_ref[i]
      target_x <- x[i]
      
      ref_a0 <- prev_table[[i]][1]
      ref_b0 <- prev_table[[i]][2]
      ref_c0 <- prev_table[[i]][3]
      ref_d0 <- prev_table[[i]][4]
      
      # contingency table of the population with asthma from a previous year
      # if ra=1, no reversibility
            a0 <- ref_b0 * (1 - ra)
            c0 <- ref_d0 * (1 - ra)
            b0 <- ref_b0 * ra
            d0 <- ref_d0 * ra
      
      # contingency table of the exposure level 
      # no exposure & no asthma: did not get asthma and did not get misdx + got asthma but misDx
            a1 <- (1 - ref_cal_inc) * ref_a0 * (1-misDx) + ref_cal_inc * ref_a0 * (1 - Dx)
      # no exposure & yes asthma: get asthma and correctly Dx + did not get asthma but misDx
            b1 <- (1 - ref_cal_inc) * ref_a0 * misDx + ref_cal_inc * ref_a0 * Dx
      # yes exposure & no asthma: got asthma but incorrectly Dx + did not get asthma and correctly Dx
            c1 <- (1 - cal_inc) * ref_c0 * (1 - misDx) + cal_inc * ref_c0 * (1 - Dx)
      # yes exposure & yes asthma: got asthma and correctly Dx
            d1 <- (1 - cal_inc) * ref_c0 * misDx + cal_inc * ref_c0 * Dx
      
      # two targets
            # objective: asthma prev OR
      a <- a0 + a1
      b <- b0 + b1
      c <- c0 + c1
      d <- d0 + d1
            tmp_OR <- a * d / (b * c)
      result <- result +
        abs(log(target_OR_no_ref[i]) - (log(d) + log(a) - log(b) - log(c)))
    }
    
        return(result %>% unlist() %>% mean())
  }
  
  fnc_value <- obj_function(log_inc_OR)
  return(fnc_value)
}


inc_correction_calculator <- function(target_inc,
                              past_target_prev,
                              past_target_OR,
                              target_OR,
                              p_risk,
                              ra=1,
                              misDx=0,
                              Dx=1,
                              risk_set,
                              log_inc_OR){
  
  beta0 <- logit(target_inc)
  
  prop_asthma <- past_target_prev # (ref_b0+ref_d0)
  prop_no_asthma <- (1-past_target_prev) # (ref_a0+ref_c0)
  
  # reconstruct contingency table for each OR
  ref_p_risk <- p_risk[1]
  sol <- prev_calibrator(past_target_prev,past_target_OR, p_risk)
  p0 <- inverse_logit(logit(past_target_prev) - sum(p_risk[-1]*sol))
  calibrated_p <- inverse_logit(logit(p0) + log(past_target_OR))
  # distribution of the risk factors for the population without asthma
  no_asthma_p_risk_dist <- (1-calibrated_p) * p_risk
  # normalize
  no_asthma_p_risk_dist <- no_asthma_p_risk_dist/sum(no_asthma_p_risk_dist)
  
  # for each OR, we need to obtain the contingency table
  
  prev_table <- c()
  
  for(i in 1:(length(past_target_OR)-1)){
    # print(i)
    tmp_p_risk <- p_risk[c(1,i+1)]
    tmp_p_risk <- tmp_p_risk/sum(tmp_p_risk)
    tmp_p <- calibrated_p[c(1,i+1)]
    # return: a b c d
    # a: no exp, no asthma
    # b: no exp, yes asthma
    # c: yes exp, no asthma
    # d: yes exp, yes asthma
    # Solve the following:
    # tmp_target_prev  = (b+d)/(a+b+c+d) 
    # tmp_p[1] = b/(a+b) 
    # tmp_p[2] = d/(c+d) 
    # tmp_target_OR  = (a*d)/(b*c) 
    
    # someone else has done it;
    # use the metafor pkg
    # Bonett, D. G. (2007).
    # Transforming odds ratios into correlations for meta-analytic research. 
    # American Psychologist, 62(3), 254–255. ⁠https://doi.org/10.1037/0003-066x.62.3.254⁠
    
    nn <- 1e10
    prev_table[[i]] <-  rev(metafor::conv.2x2(ori=past_target_OR[i+1],
                                              ni = nn,
                                              # prev of exposure
                                              n1i = ((1-tmp_p[2])*tmp_p_risk[2] +tmp_p_risk[2] * tmp_p[2])*nn,
                                              # prev of asthma
                                              n2i=  sum(tmp_p_risk * tmp_p)*nn)/nn)
  }
  
  future_prev_table <- c()
  target_prev <- (past_target_prev*ra + target_inc*(1-past_target_prev)*Dx +
                    (1-target_inc)*(1-past_target_prev)*misDx)
  tmp_sol <- prev_calibrator(target_prev,target_OR, p_risk)
  tmp_p0 <- inverse_logit(logit(target_prev) - sum(p_risk[-1]*tmp_sol))
  tmp_calibrated_p <- inverse_logit(logit(tmp_p0) + log(target_OR))
  
  for(i in 1:(length(past_target_OR)-1)){
    # print(i)
    tmp_p_risk <- p_risk[c(1,i+1)]
    tmp_p_risk <- tmp_p_risk/sum(tmp_p_risk)
    tmp_p <- tmp_calibrated_p[c(1,i+1)]
    
    nn <- 1e10
    future_prev_table[[i]] <-  rev(metafor::conv.2x2(ori=target_OR[i+1],
                                                     ni = nn,
                                                     # prev of exposure
                                                     n1i = ((1-tmp_p[2])*tmp_p_risk[2] +tmp_p_risk[2] * tmp_p[2])*nn,
                                                     # prev of asthma
                                                     n2i=  sum(tmp_p_risk * tmp_p)*nn)/nn)
  }
  
  obj_function <- function(y){
    x <- y
    
    # # x = log(OR) for incidence eqn
    
    q <- length(no_asthma_p_risk_dist)-1
    target_OR_no_ref <- target_OR[-1]
    # # calibrate the current inc to the target inc
    tmp_sol <- prev_calibrator(
        asthma_prev_target=target_inc,
        target_OR=exp(c(0,x)),
        risk_factor_prev=no_asthma_p_risk_dist
    )
    logit_p0 <- beta0 - sum(tmp_sol*no_asthma_p_risk_dist[-1])
    # logit_p0 <- beta0
    calibrated_inc <- inverse_logit(logit_p0 + c(0,x))
    
    inc_correction_term <- -sum(tmp_sol*no_asthma_p_risk_dist[-1])
    
    ref_cal_inc <- calibrated_inc[1]
    calibrated_inc_no_ref <- calibrated_inc[-1]
    
    result <- 0
    
    for(i in 1:(length(target_OR)-1)){
      cal_inc <- calibrated_inc_no_ref[i]
      target_x <- x[i]
      
      ref_a0 <- prev_table[[i]][1]
      ref_b0 <- prev_table[[i]][2]
      ref_c0 <- prev_table[[i]][3]
      ref_d0 <- prev_table[[i]][4]
      
      # contingency table of the population with asthma from a previous year
      # if ra=1, no reversibility
      a0 <- ref_b0*(1-ra)
      c0 <- ref_d0*(1-ra)
      b0 <- ref_b0*ra
      d0 <- ref_d0*ra
      
      # contingency table of the exposure level 
      # no exposure & no asthma: did not get asthma and did not get misdx + got asthma but misDx
      a1 <-  (1-ref_cal_inc)*ref_a0*(1-misDx) + ref_cal_inc * ref_a0 *(1-Dx)
      # no exposure & yes asthma: get asthma and correctly Dx + did not get asthma but misDx
      b1 <- (1-ref_cal_inc)*ref_a0*misDx + ref_cal_inc * ref_a0 * Dx
      # yes exposure & no asthma: got asthma but incorrectly Dx + did not get asthma and correctly Dx
      c1 <- (1-cal_inc)*ref_c0*(1-misDx) + cal_inc * ref_c0 * (1-Dx)
      # yes exposure & yes asthma: got asthma and correctly Dx
      d1 <-   (1-cal_inc)*ref_c0*misDx + cal_inc* ref_c0 * Dx
      
      # two targets
      #  objective: asthma prev OR
      a <- a0 + a1
      b <- b0 + b1
      c <- c0 + c1
      d <- d0 + d1
      tmp_OR <- a*d/(b*c)
      result <- result +
        abs(log(target_OR_no_ref[i]) - (log(d) + log(a) - log(b) - log(c)))
      # sum(abs(c(a,b,c,d)-future_prev_table[[i]]))
      
      # abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1)) #figure out why this works
      # abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1)) #figure out why this works
      
      # result <- result + abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1) - target_x)
    }
    
    # print(tmp_sol)
    return(c(result %>% 
             unlist() %>% 
             mean(),inc_correction_term))
  }
  
  # risk_set <- OR_generator(risk_set,inc_parameters)
  
  # n_par <-   sum(risk_set[,c(1,2)] %>%
  #                  apply(.,2,function(x){length(unique(x))-1}))
  
  fnc_value <- obj_function(log_inc_OR)
  
  return(fnc_value)
  
}