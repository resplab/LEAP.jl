# source: https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1310009608

library(tidyverse)
library(here)
library(grid)
library(gridExtra)
library(mgcv)
library(ggpubr)

max_year <- 2019
DF_STATCAN <- read_csv(here("R/public_dataset", "13100096.csv"))
chosen_province <- "British Columbia"
baseline_year <- 2000


generate_plots <- function(df, linear_model) {
    lin_data <- cbind(df, fitted=linear_model$fitted.values)
    plot_list <- c()
    counter <- 1
    for(Year in 2015:2020) {
        look2 <- df %>% 
        filter(year==Year) %>% 
        mutate(
            pred=pred*100,
            prop=prop*100,
            prop_upper=prop_upper*100,
            prop_lower=prop_lower *100
        )

        plot_list[[counter]] <-   
            ggplot(data=look2, aes(x=age,y=pred))+
            geom_line(size=2, alpha=0.7) +
            ylim(c(0,20)) +
            geom_point(
                data=look2 %>% filter(!is.na(prop)),
                aes(x=age,y=prop)
            ) +
            geom_errorbar(
                data=look2 %>% filter(!is.na(prop)),
                aes(ymin=prop_lower, ymax=prop_upper, x=age),alpha=0.8,width=3
            ) +
            geom_line(
                data=lin_data %>% filter(year==Year) %>% 
                    mutate(fitted=fitted*100),
                aes(y=fitted,x=age),col='red',size=1.5,alpha=0.4
            )+
            ggtitle(Year)
            counter <- counter+1
    }
    return(plot_list)
}


survey_data_processor <- function(chosen_province) {
    # Rename columns
    df <- DF_STATCAN %>%
        rename(
            year=REF_DATE,
            sex=Sex,
            age_group=`Age group`,
            prop=VALUE,
            characteristics=Characteristics,
            province=GEO
        )
    # Filter for the year, chosen province, characteristics, indicators, and remove age group "Total"
    df <- df %>% 
        filter(
            str_detect(province, chosen_province) &
            str_detect(characteristics, "Percent|percent") &
            !str_detect(age_group, "Total") &
            Indicators == 'Asthma' &
            year <= max_year
        )
    # Select year, characteristics, age_group, sex, and prop columns
    df <- df %>%
        select(year, characteristics, age_group, sex, prop)
    df <- df %>%
        mutate(
            sex=str_sub(sex,1,1),
            prop=prop/100,
            age_group=str_remove(age_group," years and over| years")
        )
    age_ranges <- do.call(rbind, str_split(df$age_group, ' to '))
    df$age_lower <- as.numeric(age_ranges[,1])
    df$age_upper <- as.numeric(age_ranges[,2])
    # If the upper age is the same as the lower age, set the upper age to 100
    # (this is the 65+ category)
    df <- df %>% 
        mutate(age_upper=if_else(age_lower==age_upper, 100, age_upper))
    
    # Convert 95% upper and lower intervals to columns
    df <- df %>%
        pivot_wider(names_from=characteristics, values_from=prop) %>% 
        rename(
            prop=Percent,
            prop_lower=`Low 95% confidence interval, percent`,
            prop_upper=`High 95% confidence interval, percent`
        )
  
    # Calculate the standard deviation, variance, and the alpha and beta parameters
    df <- df %>%
        mutate(
            sd=(prop_upper-prop)/qnorm(0.975),
            var=sd^2,
            sanity=prop+sd*qnorm(0.975),
            alpha=((1-prop)/var - 1/prop)*prop^2,
            beta=alpha * (1/prop-1)
        )
        
    # Remove age_group and sanity columns
    df <- df %>% select(-age_group, -sanity)
    
    # Remove "Both" sexes from the data
    df <- df %>% 
        filter(sex!="B")

    # Set the age group to the middle of the range
    df <- df %>%
        mutate(age=floor((age_lower + age_upper) / 2)
    )
  
    # CHILD STUDY DATA (restricted access)
    
    df_child <- readxl::read_xlsx(
        here("R/private_dataset/child_dataset/AJ_Cohort_Creation","CHILD_Master.xlsx"),
        1
    )
    # Column 1: Study_ID
    # Columns with "asthma" in the name
    # Column 147: Birth.Baby.s.sex
    df_child <- df_child[,c(1, which(str_detect(colnames(df_child), "asthma")), 147)]

    # Column 1: Study_ID
    # Column 3: 5.Years.Has.your.child.EVER.been.diagnosed.with.asthma.
    # Column 27: Birth.Baby.s.sex
    df_child <- df_child[, c(1,3,27)]
    with(
        df_child %>% na.omit(),
        table(`5.Years.Has.your.child.EVER.been.diagnosed.with.asthma.`, Birth.Baby.s.sex)
    )
    # Asthma    M      F
    # No      1072  1004
    # Yes      90    54

    n_asthma_M <- 90
    n_asthma_F <- 54
    n_no_asthma_M <- 1072
    n_no_asthma_F <- 1004

    144/(144+2076) # 6.49%; 5 years prev
    table(df_child$Birth.Baby.s.sex) # 1: Male; 2: Female
  
    # 90/(1072+90) = 0.07745267 (expectation value)
    sd_M <- sqrt(n_asthma_M)* qnorm(0.975)
    
    # 54/(1004+54) # 0.05104
    sd_F <-  sqrt(n_asthma_F)* qnorm(0.975)
    
    M5 <- data.frame(
        year=2015,
        sex="M",
        age=5,
        prop=n_asthma_M/(n_no_asthma_M + n_asthma_M),
        prop_upper = 1/(n_no_asthma_M + n_asthma_M)*(n_asthma_M+sd_M),
        prop_lower = 1/(n_no_asthma_M + n_asthma_M)*(n_asthma_M-sd_M)
    )
  
    F5 <- data.frame(
        year=2015,
        sex="F",
        age=5,
        prop=n_asthma_F/(n_no_asthma_F + n_asthma_F),
        prop_upper = 1/(n_no_asthma_F + n_asthma_F)*(n_asthma_F+sd_F),
        prop_lower = 1/(n_no_asthma_F + n_asthma_F)*(n_asthma_F-sd_F)
    )
    
    # do a spline stratified by sex
  
    df_female <- df %>% 
        select(year,sex,age,prop,prop_lower,prop_upper) %>% 
        # rbind(M5) %>% 
        rbind(F5) %>% 
        filter(sex=="F")
    
    spline_female <- mgcv::gam(prop~s(age, bs='cs', k=6) + year, data=df_female)
    
    tmp_data_female <- expand.grid(
        year=2015:2025,
        age=3:110,
        sex=c("F")
    )
  
    tmp_data_female$pred <- predict(spline_female, newdata=tmp_data_female)
    
    tmp_data_female <- tmp_data_female %>% 
        left_join(
            df %>% select(year,age,sex,prop,prop_upper,prop_lower) %>% rbind(F5),
            by=c("year","age","sex")
        ) %>% 
        select(-sex)
    
    linear_model_female <- stats::lm(
        pred~age+I(age^2)+I(age^3)+I(age^4)+I(age^5)+I(age^6)+I(age^7)+I(age^8)+I(age^9)+I(age^10)+I(age^11)+
            I(age^12)+year,
        data=tmp_data_female
    )

    df_male = df %>%
        select(year,sex,age,prop,prop_lower,prop_upper) %>%
        rbind(M5) %>%
        # rbind(F5) %>%
        filter(sex=="M")
  
    spline_male <- mgcv::gam(prop~s(age, bs="ps", m=1) + year, data=df_male)
  
    tmp_data_male <- expand.grid(year=2015:2025,
                                age=3:110,
                                sex=c("M"))
  
    tmp_data_male$pred <- predict(spline_male,newdata = tmp_data_male)
    
    tmp_data_male <- tmp_data_male %>% 
        left_join(
            df %>% select(year,age,sex,prop,prop_upper,prop_lower) %>% rbind(M5),
            by=c("year","age","sex")
        ) %>% 
        select(-sex)
    
    linear_model_male <- lm(
        pred~year+age+I(age^2)+I(age^3)+I(age^4)+I(age^5)+I(age^6)+I(age^7)+I(age^8)+I(age^9)+I(age^10)+I(age^11)+
            I(age^12),
        data=tmp_data_male
    )

    # Get Coefficients

    summ_female <- summary(linear_model_female)
    female_coef <- summ_female$coefficients %>% 
        as.data.frame() %>%  
        select(1,2) %>%
        mutate(sex="F")
    colnames(female_coef)[1:2] <- c("est","std")
  
    female_coef <- female_coef %>% 
        rbind(sigma=data.frame(est=summ_female$sigma, std=0, sex="F"))
    
    summ_male <- summary(linear_model_male)
    male_coef <- summ_male$coefficients %>% 
        as.data.frame() %>%  
        select(1,2) %>%
        mutate(sex="M")
    colnames(male_coef)[1:2] <- c("est","std")
    male_coef <- male_coef %>%
        rbind(sigma=data.frame(est=summ_male$sigma, std=0, sex="F"))

    CA_prev_coef <- rbind(female_coef, male_coef)

    # Create Plots

    plot_list_female <- generate_plots(df=tmp_data_female, linear_model=linear_model_female)
    plot_list_male <- generate_plots(df=tmp_data_male, linear_model=linear_model_male)

    plot_list <- c()
    for(j in 1:length(plot_list_male)) {
        plot_list[[(2*(j-1)+1)]] <- plot_list_female[[j]]+ xlab("") + ylab("") + ggtitle("") + theme(text=element_text(size=15))
        plot_list[[(2*(j-1)+2)]] <- plot_list_male[[j]] + xlab("") + ylab("") + ggtitle("") + theme(text=element_text(size=15))
    }
  
    # left: female; right: male
    grid.arrange(
        arrangeGrob(plot_list[[1]], plot_list[[2]], top="2015",nrow=1),
        # arrangeGrob(plot_list[[3]], plot_list[[4]], top="2016",nrow=1),
        # arrangeGrob(plot_list[[5]], plot_list[[6]], top="2017",nrow=1),
        # arrangeGrob(plot_list[[7]], plot_list[[8]], top="2018",nrow=1),
        # arrangeGrob(plot_list[[9]], plot_list[[10]], top="2019",nrow=1),
        arrangeGrob(plot_list[[11]], plot_list[[12]], top="2020",nrow=1),
        ncol=1,
        bottom=textGrob("Age (years)", gp=gpar(fontsize=20)),
        left=textGrob("Asthma prevalence (per 100)", rot="90", gp=gpar(fontsize=20))
    )

    return(list(
        plot_list_male=plot_list_male,
        plot_list_female=plot_list_female,
        prev_coef=CA_prev_coef
    ))
}

# BC ----------------------------------------------------------------------
tmp <- readxl::read_xlsx("private_dataset/asthma_inc_prev.xlsx",sheet=1) %>% 
    filter(age_group_desc != "<1 year") %>% 
    mutate(year=substr(fiscal_year,1,4) %>% as.numeric()) %>% 
    rename(sex=gender) %>% 
    filter(year>=baseline_year) %>% 
    rename(age_group = age_group_desc)
  
lapply(tmp$age_group,function(x){
  ceiling(mean(parse_number(str_split(x,"-")[[1]])))
})  %>% unlist() -> tmp$age

# Key assumption: set the incidence level at age 3 to the prevelance level
tmp <- tmp %>% 
    mutate(incidence = ifelse(age==3,prevalence,incidence)) %>% 
    mutate(age= ifelse(age==90,100,age))

plotter <- function(SEX, years){
  df <- tmp %>%
    filter(sex==SEX) %>%
    select(-sex)

    spline_model <- gam(prevalence~s(age,year,bs='gp',k=50) ,data=df)

    tmp_data <- expand.grid(year=years, age=3:100) %>% 
      mutate()

    tmp_data$pred <- predict(spline_model,newdata = tmp_data)


    linear_model <- lm(
        pred~age+year+I(age^2)+I(age^3)+I(age^4)+I(age^5)+I(age^6)+I(age^7)+I(age^8)+I(age^9)+I(age^10)+I(age^11)+
        I(age^12),
        data=tmp_data
    )

    lin_data <- cbind(tmp_data, fitted=linear_model$fitted.values)

    tmp_data %>%
        left_join(
            df %>% select(year,age,prevalence,prev_lower,prev_upper),
            by=c("year","age")
        ) -> tmp_data

    results <- c()
    counter <- 1
    for(Year in years){
      
      look2 <- tmp_data %>% 
        filter(year==Year) %>% 
        mutate(pred = pred*100,
               prevalence=prevalence*100,
               prev_lower = prev_lower *100,
               prev_upper = prev_upper*100)
      
      
      results[[counter]] <-
        ggplot(data=look2,aes(x=age,y=pred))+
        geom_line(size=2,alpha=0.7) +
        ylim(c(0,20))+
        geom_point(data=look2, aes(x=age,y=prevalence)) +
        geom_errorbar(data=look2,
                      aes(ymin=prev_lower,ymax=prev_upper,x=age),alpha=0.8,width=3) +
        # geom_line(data=lin_data %>% filter(year==Year),aes(y=fitted,x=age),col='red',size=1.5,alpha=0.4)+
        ggtitle(Year)
      counter <- counter+1
    }
    results
}


BC_survey_results <- survey_data_processor("British Columbia")
CA_survey_results <- survey_data_processor("Canada")

# obtain multiplier for each age; averaged across years
tmp <- c()
for(yr in 1:length(CA_plots[[1]])){
    tmp[[yr]] <- cbind(
        BC_survey_results$plot_list_male[[yr]]$data[,3]/CA_survey_results$plot_list_male[[yr]]$data[,3],
        BC_survey_results$plot_list_female[[yr]]$data[,3]/CA_survey_results$plot_list_female[[yr]]$data[,3])
}

multiplier <- 1/apply(simplify2array(tmp), 1:2, mean)


# BC asthma prev and inc equation -----------------------------------------
master_BC_asthma <- readxl::read_xlsx(here("R/private_dataset/asthma_inc_prev.xlsx"), sheet=1) %>% 
  filter(age_group_desc != "<1 year") %>% 
  mutate(year=substr(fiscal_year,1,4) %>% as.numeric()) %>% 
  rename(sex=gender) %>% 
  filter(year>=1996) %>% 
  rename(age_group = age_group_desc) %>% 
  mutate(prev_sd = sqrt(prevalence_numerator)*qnorm(0.975),
         prev_upper = (prevalence_numerator+prev_sd)/pop,
         prev_lower = (prevalence_numerator-prev_sd)/pop) 

lapply(master_BC_asthma$age_group,function(x){
  ceiling(mean(parse_number(str_split(x,"-")[[1]])))
})  %>% 
  unlist() -> master_BC_asthma$age

# Key assumption: asthma starts at age 3
# Set the inc = prev at age 3
master_BC_asthma <- master_BC_asthma %>% 
  mutate(incidence = ifelse(age==3,prevalence,incidence),
         incidence_numerator = ifelse(age==3,prevalence_numerator,incidence_numerator)) 

# Male
plotter_sex_prev <- function(SEX,df_master,pp,years,ylim_lower=0,ylim_upper=22){
  df <- df_master %>% 
    filter(sex==SEX)
  
  tmp_data <- c()
  
  unique_years <- as.integer(years) %>% 
    unique()
    
  for(YEAR in unique_years){
    
    spline_model <- gam(y~s(age,bs='ps',m=1),data=df %>% 
                          filter(year==YEAR))
    
    tmptmp_data <- expand.grid(year=YEAR, age=seq(3,95,by=0.01)) %>% 
      as.data.frame()
    
    tmptmptmp <- predict(spline_model,
                         tmptmp_data %>% 
                           select(-year),se.fit=T)
    tmptmp_data$pred <- tmptmptmp$fit
    tmptmp_data$se <- tmptmptmp$se.fit
    tmp_data <- rbind(tmp_data,tmptmp_data)
  }
  
 
  # tmp_data %>% filter(age%in% c(3:95)) %>% 
  #   left_join(df,by=c('year',"age")) %>% 
  #   filter(!is.na(prevalence))
  #
  plotter <- function(years,df,tmp_data){
    
    # tmp_data$actual <- df$prevalence
    tmptmp_data <- tmp_data %>% 
      mutate(year=year/1000,
             age=age/100)

    # tmptmp_dat
    # pp <- 5
    tmp_age <- model.matrix(~poly(tmptmp_data$age,degree=pp,raw=T))[,-c(1:2)]
    colnames(tmp_age) <- paste0("age",c(2:pp))
    tmp_year <-  model.matrix(~poly(tmptmp_data$year,degree=pp,raw=T))[,-c(1:2)]
    colnames(tmp_year) <- paste0("year",c(2:pp))
    tmp_age_year <- tmp_age*tmp_year
    colnames(tmp_age_year) <- paste0("ageyear",c(2:pp))
    
    
    linear_model <- lm(pred~.,data=cbind(tmptmp_data,tmp_age,tmp_year,tmp_age_year))
    lin_data <- cbind(tmp_data,fitted=linear_model$fitted.values)
    
    tmp_data %>%
      left_join(df %>% select(year,age,y,y_lower,y_upper),
                by=c("year","age")) -> tmp_data
    
    results <- c()
    counter <- 1
    for(Year in (as.integer(years) %>% unique()) ){
      
      look2 <- tmp_data %>% 
        # mutate(year=year*1000,
               # age=age*100) %>% 
        filter(year==Year) %>% 
        mutate(pred = pred*100,
               se = se*100,
               y=y*100,
               y_lower = y_lower *100,
               y_upper = y_upper*100) %>% 
        mutate(pred_lower = pred - qnorm(0.975)*se,
               pred_upper = pred + qnorm(0.975)*se)
      
      
      results[[counter]] <-
        ggplot(data=look2,aes(x=age,y=pred))+
        geom_ribbon(aes(ymin=pred-qnorm(0.975)*se,ymax=pred+qnorm(0.975)*se),fill='grey70')+
        geom_line(size=2,alpha=0.7) +
        ylim(c(ylim_lower,ylim_upper))+
        geom_point(data=look2, aes(x=age,y=y)) +
        geom_errorbar(data=look2,
                      aes(ymin=y_lower,ymax=y_upper,x=age),alpha=0.8,width=3) +
        # geom_line(data=lin_data %>% filter(year==Year),aes(y=fitted*100,x=age),col='red',size=1.5,alpha=0.4)+
        # facet_grid(sex~.)+
        ggtitle(Year)
      counter <- counter+1
    }
    return(list(data=tmp_data,plots=results,model=linear_model))
    
  }

  plotter(years,df=df,tmp_data)
}

BC_asthma_prev <- master_BC_asthma %>% 
  mutate(age= ifelse(age==90,95,age)) %>% 
  select(year,sex,age,pop,prevalence,prev_sd,prev_upper,prev_lower) %>% 
  mutate(sex=as.numeric(sex=="M")) %>% 
  rename(y=prevalence,
         y_lower=prev_lower,
         y_upper=prev_upper)

# male
tmp1 <- plotter_sex_prev(1,df_master=BC_asthma_prev,pp=10,years = seq(1999,2019,by=1))
# female
tmp0 <- plotter_sex_prev(0,df_master=BC_asthma_prev,pp=10,years = seq(1999,2019,by=1))

indices <- c(2,10,15,21)

example_prev <- c(tmp1$plots[c(2,10,15,21)],tmp0$plots[c(2,10,15,21)]) %>% 
  lapply(.,function(x){
    x + 
      xlab("") +
      ylab("") +
      theme(text=element_text(size=20),
            plot.title = element_text(hjust = 0.5))+
      xlim(c(3,100))
  })
example_prev[[1]] <- example_prev[[1]] +
  ylab("Male")

example_prev[5:8] <- lapply(example_prev[5:8],
                            function(x){
                              x +
                                ggtitle("")
                            })
example_prev[[5]] <- example_prev[[5]] +
  ylab("Female")

y.grob <- textGrob("Asthma prevalence per 100 persons", 
                   gp=gpar(fontface="bold", fontsize=20), rot=90)

x.grob <- textGrob("Age (year)", 
                   gp=gpar(fontface="bold", fontsize=20))

plt <- cowplot::plot_grid(plotlist = example_prev,ncol=4)

grid.arrange(arrangeGrob(plt,left=y.grob,bottom=x.grob))

interpolated_prev <- rbind(tmp1$data %>% mutate(sex=1),
                           tmp0$data %>% mutate(sex=0))%>% 
  filter(age %in% c(3:95))

# comparison
BC_survey <- rbind(BC_survey_results$results_male[[1]]$data %>%
  filter(!is.na(prop)) %>%
  rename(prevalence=prop,
         prev_upper = prop_upper,
         prev_lower=prop_lower) %>%
    mutate(sex=1) %>% 
  select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
  mutate(DB='BC-Survey'),
  BC_survey_results$results_male[[5]]$data %>%
  filter(!is.na(prop)) %>%
  rename(prevalence=prop,
         prev_upper = prop_upper,
         prev_lower=prop_lower) %>%
  mutate(sex=1) %>% 
  select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
  mutate(DB='BC-Survey'),
  BC_survey_results$results_female[[1]]$data %>%
    filter(!is.na(prop)) %>%
    rename(prevalence=prop,
           prev_upper = prop_upper,
           prev_lower=prop_lower) %>%
    mutate(sex=0) %>% 
    select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
    mutate(DB='BC-Survey'),
  BC_survey_results$results_female[[5]]$data %>%
    filter(!is.na(prop)) %>%
    rename(prevalence=prop,
           prev_upper = prop_upper,
           prev_lower=prop_lower) %>%
    mutate(sex=0) %>% 
    select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
    mutate(DB='BC-Survey'))

CA_survey <- rbind(CA_survey_results$results_male[[1]]$data %>%
                     filter(!is.na(prop)) %>%
                     rename(prevalence=prop,
                            prev_upper = prop_upper,
                            prev_lower=prop_lower) %>%
                     mutate(sex=1) %>% 
                     select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
                     mutate(DB='CA-Survey'),
                   CA_survey_results$results_male[[5]]$data %>%
                     filter(!is.na(prop)) %>%
                     rename(prevalence=prop,
                            prev_upper = prop_upper,
                            prev_lower=prop_lower) %>%
                     mutate(sex=1) %>% 
                     select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
                     mutate(DB='CA-Survey'),
                   CA_survey_results$results_female[[1]]$data %>%
                     filter(!is.na(prop)) %>%
                     rename(prevalence=prop,
                            prev_upper = prop_upper,
                            prev_lower=prop_lower) %>%
                     mutate(sex=0) %>% 
                     select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
                     mutate(DB='CA-Survey'),
                   CA_survey_results$results_female[[5]]$data %>%
                     filter(!is.na(prop)) %>%
                     rename(prevalence=prop,
                            prev_upper = prop_upper,
                            prev_lower=prop_lower) %>%
                     mutate(sex=0) %>% 
                     select(year,sex,age,pred,prevalence,prev_lower,prev_upper) %>%
                     mutate(DB='CA-Survey'))

admin_2015_2019 <- interpolated_prev %>%
  filter(year %in% c(2015,2019)) %>% 
  filter(!is.na(y)) %>%
  mutate(DB="BC-Admin") %>% 
  rename(prevalence=y,
         prev_upper = y_upper,
         prev_lower=y_lower) %>%
  mutate(prevalence=prevalence*100,
         prev_upper = prev_upper*100,
         prev_lower= prev_lower*100) %>% 
  select(year,sex,age,pred,prevalence,prev_lower,prev_upper,DB)

df_survey_admin <- rbind(admin_2015_2019,BC_survey,CA_survey)

ggplot(df_survey_admin %>% 
         mutate(sex=factor(sex,c(0,1),c("Female","Male"))),aes(x=age,y=prevalence,color=DB)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=prev_lower,ymax=prev_upper,x=age),alpha=0.8,width=3) +
  geom_line()+
  ylab("Prevalence of asthma (per 100)") +
  xlab("Age (years)") +
  facet_grid(sex~year) + 
  theme(legend.position="top",
        legend.title = element_blank(),
        text = element_text(size=20))

############################ 
# conclusion: difference between the admin and survey data
# => use admin

BC_asthma_inc <- master_BC_asthma %>% 
  mutate(age= ifelse(age==90,95,age)) %>% 
  mutate(inc_upper = incidence+qnorm(0.975)*sqrt(incidence_numerator/pop^2),
         inc_lower =incidence-qnorm(0.975)*sqrt(incidence_numerator/pop^2)) %>% 
  select(year,sex,age,pop,incidence,inc_lower,inc_upper) %>% 
  mutate(sex=as.numeric(sex=="M")) %>% 
  rename(y=incidence,
         y_lower=inc_lower,
         y_upper=inc_upper)


plotter_sex_inc <- function(SEX,df_master,pp,years,ylim_lower=0,ylim_upper=22){
  df <- df_master %>% 
    filter(sex==SEX)
  
  tmp_data <- c()
  
  unique_years <- as.integer(years) %>% 
    unique()
  
  for(YEAR in unique_years){
    
    tmp_df <- df %>% 
      filter(year==YEAR) %>% 
      arrange(year,sex,age)
    
    spline_model <- gam(y~s(age,bs='ps'),data=tmp_df,family = gaussian(link='log'))
    
    tmptmp_data <- expand.grid(year=YEAR, age=seq(3,95,by=0.01)) %>% 
      as.data.frame()
    
    tmptmptmp <- predict(spline_model,
                         tmptmp_data %>% 
                           select(-year),se.fit=T,type='response')
    tmptmp_data$pred <- tmptmptmp$fit
    tmptmp_data$se <- tmptmptmp$se.fit
    tmp_data <- rbind(tmp_data,tmptmp_data)
  }
  
  
  # tmp_data %>% filter(age%in% c(3:95)) %>% 
  #   left_join(df,by=c('year',"age")) %>% 
  #   filter(!is.na(prevalence))
  #
  plotter <- function(years,df,tmp_data){
    
    # tmp_data$actual <- df$prevalence
    tmptmp_data <- tmp_data %>% 
      mutate(year=year/1000,
             age=age/100)
    
    # tmptmp_dat
    # pp <- 5
    tmp_age <- model.matrix(~poly(tmptmp_data$age,degree=pp,raw=T))[,-c(1:2)]
    colnames(tmp_age) <- paste0("age",c(2:pp))
    tmp_year <-  model.matrix(~poly(tmptmp_data$year,degree=pp,raw=T))[,-c(1:2)]
    colnames(tmp_year) <- paste0("year",c(2:pp))
    tmp_age_year <- tmp_age*tmp_year
    colnames(tmp_age_year) <- paste0("ageyear",c(2:pp))
    
    
    linear_model <- lm(pred~.,data=cbind(tmptmp_data,tmp_age,tmp_year,tmp_age_year))
    lin_data <- cbind(tmp_data,fitted=linear_model$fitted.values)
    
    tmp_data %>%
      left_join(df %>% select(year,age,y,y_lower,y_upper),
                by=c("year","age")) -> tmp_data
    
    results <- c()
    counter <- 1
    for(Year in (as.integer(years) %>% unique()) ){
      
      look2 <- tmp_data %>% 
        # mutate(year=year*1000,
        # age=age*100) %>% 
        filter(year==Year) %>% 
        mutate(pred = pred*100,
               se = se*100,
               y=y*100,
               y_lower = y_lower *100,
               y_upper = y_upper*100) %>% 
        mutate(pred_lower = pred - qnorm(0.975)*se,
               pred_upper = pred + qnorm(0.975)*se)
      
      
      results[[counter]] <-
        ggplot(data=look2,aes(x=age,y=pred))+
        geom_ribbon(aes(ymin=pred-qnorm(0.975)*se,ymax=pred+qnorm(0.975)*se),fill='grey70')+
        geom_line(size=2,alpha=0.7) +
        ylim(c(ylim_lower,ylim_upper))+
        geom_point(data=look2, aes(x=age,y=y)) +
        geom_errorbar(data=look2,
                      aes(ymin=y_lower,ymax=y_upper,x=age),alpha=0.8,width=3) +
        # geom_line(data=lin_data %>% filter(year==Year),aes(y=fitted*100,x=age),col='red',size=1.5,alpha=0.4)+
        # facet_grid(sex~.)+
        ggtitle(Year)
      counter <- counter+1
    }
    return(list(data=tmp_data,plots=results,model=linear_model))
    
  }
  
  plotter(years,df=df,tmp_data)
}


inc_tmp1 <- plotter_sex_inc(1,df_master=BC_asthma_inc,pp=15,years = seq(1999,2019,by=1),ylim_upper=12)
inc_tmp0 <- plotter_sex_inc(0,df_master=BC_asthma_inc,pp=15,years = seq(1999,2019,by=1),ylim_upper=8)


indices <- c(2,10,15,21)

example_inc <- c(inc_tmp1$plots[c(2,10,15,21)],inc_tmp0$plots[c(2,10,15,21)]) %>% 
  lapply(.,function(x){
    x + 
      xlab("") +
      ylab("") +
      theme(text=element_text(size=20),
            plot.title = element_text(hjust = 0.5))+
      xlim(c(3,100))
  })
example_inc[[1]] <- example_inc[[1]] +
  ylab("Male")

example_inc[5:8] <- lapply(example_inc[5:8],
                            function(x){
                              x +
                                ggtitle("")
                            })
example_inc[[5]] <- example_inc[[5]] +
  ylab("Female")

y.grob <- textGrob("Asthma incidence per 100 persons", 
                   gp=gpar(fontface="bold", fontsize=20), rot=90)

x.grob <- textGrob("Age (year)", 
                   gp=gpar(fontface="bold", fontsize=20))

plt <- cowplot::plot_grid(plotlist = example_inc,ncol=4)

grid.arrange(arrangeGrob(plt,left=y.grob,bottom=x.grob))


interpolated_inc <- rbind(inc_tmp1$data %>% mutate(sex=1),
                          inc_tmp0$data %>% mutate(sex=0)) %>% 
  filter(age %in% c(3:110))

# set prev = inc at age 3
look <- data.frame(inc3 = interpolated_inc %>% 
  filter(age==3) %>% select(pred) %>% unlist() %>% as.vector(),
  prev3 = interpolated_prev %>% filter(age==3) %>% 
    select(pred) %>% unlist() %>% as.vector()) %>% 
  mutate(diff = inc3-prev3)


summary(look$diff) # minimal

interpolated_inc$pred[which(interpolated_inc$age==3)] <- interpolated_prev$pred[which(interpolated_prev$age==3)]

# write_csv(interpolated_inc,"BC_asthma_inc_interpolated.csv")
# write_csv(interpolated_prev,"BC_asthma_prev_interpolated.csv")

# CA asthma prev
tmp_multiplier <- cbind(age=3:110,multiplier=multiplier)
# 1 is male
colnames(tmp_multiplier)[2:3] <- c("Male","Female")
tmp_multiplier %>% 
  as.data.frame() %>% 
  pivot_longer(-1,names_to="sex",values_to="multiplier") %>% 
  mutate(sex=as.numeric(sex=="Male")) -> tmp_multiplier

interpolated_prev %>% 
  left_join(tmp_multiplier,by=c('age','sex')) %>% 
  mutate(pred= pred * multiplier,
         se = se * multiplier,
         y = y * multiplier,
         y_lower = y * multiplier,
         y_upper = y_upper * multiplier) -> interpolated_prev_CA

interpolated_inc %>% 
  left_join(tmp_multiplier,by=c('age','sex')) %>% 
  mutate(pred= pred * multiplier,
         se = se * multiplier,
         y = y * multiplier,
         y_lower = y * multiplier,
         y_upper = y_upper * multiplier) -> interpolated_inc_CA

master_prev <- rbind(interpolated_prev %>% mutate(province="BC"),
                     interpolated_prev_CA %>%  select(-multiplier) %>% 
                       mutate(province="CA")) %>% 
  select(year,age,pred,sex,province) %>%
  mutate(sex = ifelse(sex==1,"M","F")) %>% 
  pivot_wider(names_from=sex,values_from=pred) %>% 
  select(year,age,`F`,M,province)


master_inc <- rbind(interpolated_inc %>% mutate(province="BC"),
                    interpolated_inc_CA %>% select(-multiplier) %>% mutate(province="CA")) %>% 
  select(year,age,pred,sex,province) %>%
  mutate(sex = ifelse(sex==1,"M","F")) %>% 
  pivot_wider(names_from=sex,values_from=pred) %>% 
  select(year,age,`F`,M,province)


write_csv(master_prev,"master_asthma_prev_interpolated.csv")
write_csv(master_inc,"master_asthma_inc_interpolated.csv")

# check

master_prev %>% filter(year==baseline_year & province=="BC")
master_inc %>% filter(year==baseline_year & province=="BC")

master_prev %>% filter(year==baseline_year & province=="CA")
master_inc %>% filter(year==baseline_year & province=="CA")