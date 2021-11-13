## Libraries we need
library(socialmixr)
library(pracma)
library(deSolve)
library(ggplot2)
library(tidyverse)
library(shiny)
library(getTBinR)
library(patchwork)

source("code/auxiliary_funcs.R")
source("code/sir_functions.R")

age_distribution <- read_csv("parameters/age_distribution.csv")
## POLYMOD contact matrix
data("polymod")
N_props1 <- age_distribution$age_dist
age_dat <- data.frame(lower.age.limit=seq(0,80,by=10),population=N_props1)
beta_scales <- age_distribution$infectivity
alphas1 <- age_distribution$susceptibility
polymod$contacts <- polymod$contacts %>% group_by(part_id) %>% sample_n(5,replace=TRUE)
polymod_c <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)

## Solve base model for comparison
N_props <- age_distribution %>% pull(age_dist)
N_age_classes <- length(N_props1)
N_immunity_classes <- 9

epi_curve <- read_csv("parameters/full_epi_curve.csv")
epi_report <- read_csv("parameters/epi_curve_2.csv")

## Number of people in each age group and age class
print("Loaded")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    output$main_plot <- renderPlot({
        parameters <- reactiveValuesToList(input)
        ts <- seq(0,parameters$tmax,by=1) ## Vector of times to solve model over. 
        
        R0 <- parameters$R0
        gamma <- parameters$gamma
        alpha <- parameters$alpha
        symptomatic_proportion <- parameters$symptomatic_proportion
        reporting_rate <- parameters$reporting_rate
        N_tot <- parameters$N_tot
        
        beta_scales[1] <- parameters$age_1_infect
        alphas1[1] <- parameters$age_1_sus
        beta_scales[2] <- parameters$age_2_infect
        alphas1[2] <- parameters$age_2_sus
        beta_scales[3] <- parameters$age_3_infect
        alphas1[3] <- parameters$age_3_sus
        beta_scales[4] <- parameters$age_4_infect
        alphas1[4] <- parameters$age_4_sus
        beta_scales[5] <- parameters$age_5_infect
        alphas1[5] <- parameters$age_5_sus
        beta_scales[6] <- parameters$age_6_infect
        alphas1[6] <- parameters$age_6_sus
        beta_scales[7] <- parameters$age_7_infect
        alphas1[7] <- parameters$age_7_sus
        beta_scales[8] <- parameters$age_8_infect
        alphas1[8] <- parameters$age_8_sus
        beta_scales[9] <- parameters$age_9_infect
        alphas1[9] <- parameters$age_9_sus
        
        print(beta_scales)
        print(alphas1)
        
        I0 <- parameters$I0 ## Proportion of population initially infected, the seed
        
        N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)
        
        N_props_long <- N_props
        N <- matrix(N_props_long*N_tot,ncol=1,nrow=N_age_classes)
        C_use <- setup_C(C, N, beta_scales)
        beta_par <- get_beta_vector(C,age_dat$population,gamma, R0,beta_scales)
        get_R0_vector(C,age_dat$population,gamma, R0,beta_scales)
        
        y_base <- epi_ode_size(C_use, beta_par, gamma, Ta=alpha, N, ts=ts,
                               alphas=alphas1, age_seed=4,immunity_seed=1,return_full=TRUE,seed_size=I0)
        y_end <- epi_ode_size(C_use, beta_par, gamma, Ta=alpha, N, ts=ts,
                               alphas=alphas1, age_seed=4,immunity_seed=1,return_full=FALSE,seed_size=I0)
        print(y_end)
        incidence <- y_base[,which(colnames(y_base)=="new_infs")]
        incidence <- as.data.frame(apply(incidence, 2, function(x) c(0,diff(x))))
        print(sum(incidence))
        incidence$t <- ts
    
        colnames(incidence) <- c(age_distribution$age_group,"t")
        
        incidence_melted <- reshape2::melt(incidence,id.vars="t")
        colnames(incidence_melted) <- c("t","age_group","inc")
        incidence_melted$ver <- "Infections"
        incidence_melted$report_prob <- 1
        
        symptom_onsets <- y_base[,which(colnames(y_base)=="new_symptoms")]
        symptom_onsets <- as.data.frame(apply(symptom_onsets, 2, function(x) c(0,diff(x))))
        symptom_onsets$t <- ts
        colnames(symptom_onsets) <- c(age_distribution$age_group,"t")
        symptom_onsets_melted <- reshape2::melt(symptom_onsets,id.vars="t")
        colnames(symptom_onsets_melted) <- c("t","age_group","inc")
        symptom_onsets_melted$ver <- "Symptom onsets"
        symptom_onsets_melted$report_prob <- reporting_rate*symptomatic_proportion
        
        all_cases <- bind_rows(incidence_melted,symptom_onsets_melted) %>% mutate(inc_report=inc*report_prob)
        all_cases$date <- as.Date(all_cases$t,origin="2022-07-01")
        
        p1 <- ggplot(all_cases %>% filter(ver=="Infections") %>% 
                         group_by(date,ver) %>% summarize(y=sum(inc_report)) %>% rename(Type=ver)) + 
            #geom_bar(data=epi_report,aes(x=report_date,y=n,fill=`Age group`),stat="identity") +
            geom_line(data=epi_curve,aes(x=onset_date,y=n,col="Simulation truth"),size=1.25) +
            geom_line(aes(x=date,y=y,col=`Type`),size=1) +
            scale_color_manual(name="Model prediction", values=c("Infections"="red","Symptom onsets"="orange","Simulation true infections"="black")) +
            scale_fill_who(palette="main") +
            theme_bw() +
            theme(legend.position="bottom",
                  axis.title=element_text(size=14),
                  axis.text=element_text(size=12),
                  legend.title=element_text(size=14),
                  legend.text=element_text(size=12)) +
            xlab("Date") +
            ylab("Infections") +
            guides(fill=guide_legend(nrow=3),
                   color=guide_legend(nrow=3))
        
        
        p2 <- ggplot(all_cases %>% filter(ver=="Symptom onsets") %>% 
                         group_by(date,ver) %>% summarize(y=sum(inc_report)) %>% rename(Type=ver)) + 
            geom_bar(data=epi_report,aes(x=report_date,y=n,fill=`Age group`),stat="identity") +
            #geom_line(data=epi_curve,aes(x=onset_date,y=n,col="Simulation truth"),size=1) +
            geom_line(aes(x=date,y=y,col=`Type`),size=1.25) +
            scale_color_manual(name="Model prediction", values=c("Infections"="red","Symptom onsets"="orange","Simulation true infections"="black")) +
            scale_fill_who(palette="main") +
            theme_bw() +
            theme(legend.position="bottom",
                  axis.title=element_text(size=14),
                  axis.text=element_text(size=12),
                  legend.title=element_text(size=14),
                  legend.text=element_text(size=12)) +
            xlab("Date") +
            ylab("Reported cases") +
            guides(fill=guide_legend(nrow=3),
                   color="none")
        
        age_dist_model <- all_cases %>% filter(ver=="Infections") %>% 
            group_by(age_group) %>% 
            summarize(y=sum(inc_report)) %>%
            ungroup() %>%
            mutate(prop=y/sum(y))
        
        epi_report$`Age group` <- ordered(epi_report$`Age group`,levels=age_distribution$age_group)
        age_dist_obs <- epi_report %>% 
            group_by(`Age group`) %>% 
            summarize(y=sum(n)) %>%
            ungroup() %>%
            mutate(prop=y/sum(y))
        
        p3 <- ggplot(age_dist_model) +
            geom_bar(aes(x=factor(age_group),y=prop,fill="Truth"),stat="identity")+
            geom_line(data=age_dist_obs,aes(x=factor(`Age group`),y=prop,group=1,col="Observed"),size=1)+
            scale_fill_manual(name="",values=c("Truth"="grey40")) +
            scale_color_manual(name="",values=c("Observed"="blue")) +
            theme_bw() + 
            theme(legend.position="bottom",
                  axis.title=element_text(size=14),
                  axis.text=element_text(size=12),
                  legend.title=element_text(size=14),
                  legend.text=element_text(size=12)) +
            xlab("Age group") +
            ylab("Proportion of cases")
        
        p_top <- (p1 + p2) + plot_layout(nrow=1)
        p_bot <- (p3 + plot_spacer()) + plot_layout(nrow=1)
        (p_top / p_bot) + plot_layout(heights=c(1,1),nrow=2)
        
    },height=800,width=1000)
    
    
})
