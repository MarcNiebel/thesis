#Created by Marc Niebel April 2021
#Purpose is to carry out survival analysis on acute metadata
#taking into consideration re-infection at univariable level

#Libraries needed
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)

#Variables and time_1,time2 to event
variable_data <- read.csv("Output/survival analysis_uni_sc/survival_df.csv")

#Cleaning up what  is not needed
variable_data_cleaned <- variable_data  %>%
    select(-c(spont_clearance,multiple_genotypes,acute_infection))

#Making age dichotmous
variable_data_cleaned <- variable_data_cleaned %>% mutate(Age=case_when(Age <46 ~ "<=45",
                                                                       TRUE ~ ">45"))
variable_data_cleaned$Age <- as.factor(variable_data_cleaned$Age)

#Change of reference
variable_data_cleaned$Gender <- relevel(variable_data_cleaned$Gender,ref = "Male")
variable_data_cleaned$hiv_pos <- relevel(variable_data_cleaned$hiv_pos,ref = "Yes")
variable_data_cleaned$core_hbv_Ab <-relevel(variable_data_cleaned$core_hbv_Ab,ref="Negative")
variable_data_cleaned$Cocaine_use <- relevel(variable_data_cleaned$Cocaine_use,ref = "No")
variable_data_cleaned$Heroin_use <- relevel(variable_data_cleaned$Heroin_use,ref="No")
variable_data_cleaned$Methamphetamine_use <- relevel(variable_data_cleaned$Methamphetamine_use,ref="No")
variable_data_cleaned$PWID <- relevel(variable_data_cleaned$PWID,ref = "No")
variable_data_cleaned$MSM <- relevel(variable_data_cleaned$MSM,ref = "Yes")

#Combining IL28B levels
variable_data_cleaned <- variable_data_cleaned %>% 
    mutate(IL28B_genotype_binary=case_when(IL28B_genotype == "CT"| 
                                               IL28B_genotype == "TT" ~"non-CC",
                                           IL28B_genotype =="CC" ~ "CC",
                                           TRUE ~ NA_character_)) %>%
    select(-IL28B_genotype)
variable_data_cleaned$IL28B_genotype_binary <- relevel(factor(variable_data_cleaned$IL28B_genotype_binary),ref = "non-CC")

#HIV positive only
variable_data_cleaned$Group_CD4_count <-relevel(factor(variable_data_cleaned$Group_CD4_count),ref=">=200")
variable_data_cleaned$hiv_treatment_at_hcv_diag. <- relevel(variable_data_cleaned$hiv_treatment_at_hcv_diag.,ref="No")


#Modification of dataframe to adopt a PWP gap-time model
#Getting the difference
variable_data_cleaned_diff <- variable_data_cleaned %>%
    mutate(diff=time_2-time_1) %>%
    select(-time_2)

#Assigning time_1 to 0
variable_data_cleaned_diff <- variable_data_cleaned_diff %>%
    mutate(time_1="0") %>%
    mutate_at(vars(time_1),~as.numeric(as.character(.)))

#Rename time_2
variable_data_cleaned_diff <- rename(variable_data_cleaned_diff,time_2=diff)

#Enumerate each patient
variable_data_cleaned_diff_enum <- variable_data_cleaned_diff %>%
    group_by(record_id) %>%
    mutate(enum=row_number())

#Cox proportional hazards regression####
#List of variables of interest(dichotomous ones)
dependent_variables <- c("hiv_pos","Age","Gender","PWID","core_hbv_Ab","hiv_treatment_at_hcv_diag.",
                         "peak_bilirubin_binary","Binary_peakALT","Group_CD4_count",
                         "Cocaine_use","Heroin_use","Methamphetamine_use","IL28B_genotype_binary")

#Creating formulas for each variable in a list
cox_formulas <- sapply(dependent_variables,
                       function(var) as.formula(paste('Surv(time_1,time_2,Event)~',var,'+cluster(record_id)+strata(enum)')))

#Creating cox models for each variable
cox_models <- lapply(cox_formulas,function(model){coxph(model,data=variable_data_cleaned_diff_enum)})

#Extracting the relevant data which is hazard ratio, 95% CI and p-value
cox_results <- lapply(cox_models,
                      function(rel_results){
                          rel_results<-summary(rel_results)
                          p.value <- signif(rel_results$wald["pvalue"],digits = 2)
                          hazard_ratio<-signif(rel_results$coef[2],digits = 2)
                          confint.lower <- signif(rel_results$conf.int[,"lower .95"],2)
                          confint.higher <- signif(rel_results$conf.int[,"upper .95"],2)
                          update_HR <- paste0(hazard_ratio, " (",
                                              confint.lower, "-", confint.higher,")")
                          combined_result<- c(update_HR,p.value)
                          names(combined_result) <- c("HR(95% CI)","p-value")
                          return(combined_result)
                      })
#Combining results from a list into a dataframe
cox_results_df <-do.call(rbind,lapply(cox_results,function(x)as.data.frame(t(x))))
write.csv(cox_results_df,"Output/survival analysis_uni_sc/cox_dichotomous_variables.csv")

#Cox results for multi-level variable i.e. infected genotype
#Infected genotype
cox_genotype <-coxph(Surv(time_1,time_2,Event)~combined_genotype+cluster(record_id)+strata(enum),data=variable_data_cleaned_diff_enum)
cox_genotype_output <- summary(cox_genotype)
sink("Output/survival analysis_uni_sc/cox_infected_genotype.txt")
print(cox_genotype_output)
sink(file=NULL)

#Proportional hazards assumption test####
#Using Schoenfeld residuals which should be independent of time
sch_residuals_dichotomous <- lapply(cox_models, function(x)cox.zph(x,global = FALSE))
sch_residuals_genotype <- cox.zph(cox_genotype)


#Graphical diagnostics using ggcoxzph function which includes schoenfeld results
graphical_dichotomous <- lapply(sch_residuals_dichotomous,function(x)ggcoxzph(x))
graphical_genotype <- ggcoxzph(sch_residuals_genotype)

#Outputs of both
pdf("Output/survival analysis_uni_sc/dichotomous_sch_residuals.pdf")
print(graphical_dichotomous)
dev.off()

pdf("Output/survival analysis_uni_sc/genotype_sch_residuals.pdf")
print(graphical_genotype)
dev.off()

#Incidence rate of spontaneous clearance####
#Includes both primary and reinfection
events_sc <- sum(variable_data_cleaned_diff_enum$Event)
days_at_risk <- sum(variable_data_cleaned_diff_enum$time_2)
rate_time <- (events_sc)/(days_at_risk)#cases per person-day
rate_time_person_years <- rate_time*365*100

#Primary spontaneous clearance
primary_spontaneous_clearance <- variable_data_cleaned_diff_enum %>%
    filter(enum=="1")
events_primary_sc <- sum(primary_spontaneous_clearance$Event)
days_at_risk_primary_sc <- sum(primary_spontaneous_clearance$time_2)
rate_time_primary_sc <- (events_primary_sc)/(days_at_risk_primary_sc)#cases per person-day
rate_time_person_years_primary_sc <- rate_time_primary_sc*365*100

#Reinfection spontaneous clearance
reinfection_spontaneous_clearance <- variable_data_cleaned_diff_enum %>%
    filter(enum=="2")
events_reinfection_sc <- sum(reinfection_spontaneous_clearance$Event)
days_at_risk_reinfection_sc <- sum(reinfection_spontaneous_clearance$time_2)
rate_time_reinfection_sc <- (events_reinfection_sc)/(days_at_risk_reinfection_sc)#cases per person-day
rate_time_person_years_reinfection_sc <- rate_time_reinfection_sc*365*100

#Days to spontaneous clearance in primary infection####
days_sc <- variable_data_cleaned_diff_enum %>% filter(Event == "1" & enum =="1")
#Median time to sc in days
median(days_sc$time_2)
#Quantiles in days
quantile(days_sc$time_2)
#Taking 3 months being 91 days, 6 months being 183 and 12 months being 365
time_to_sc <- as.data.frame(days_sc$time_2)
time_to_sc_cal <- time_to_sc %>%
    mutate(time_grouping=case_when(`days_sc$time_2` < 92 ~ "3 months",
                                   `days_sc$time_2` >91  & `days_sc$time_2` < 184 ~ "6 months",
                                   `days_sc$time_2` > 183  & `days_sc$time_2`  <366  ~ "12 months",
                                   TRUE ~ ">1 year"))
#This gives an overview of spontaneous clearance time
time_to_sc_cal %>%
    group_by(time_grouping) %>%
    summarise(n=n())

#Kaplan-Meier plots####
#Only of relevance for single event analysis
KM_data_primary_event <- variable_data_cleaned_diff_enum %>%
    filter(enum=="1")

#Make a list of of formulas
formula_age <- survfit(Surv(time_2,Event)~Age, data=KM_data_primary_event)
formula_gender <- survfit(Surv(time_2,Event)~Gender, data=KM_data_primary_event)
formula_hiv <- survfit(Surv(time_2,Event)~hiv_pos, data=KM_data_primary_event)
formula_pwid <- survfit(Surv(time_2,Event)~PWID, data=KM_data_primary_event)
formula_hbv <- survfit(Surv(time_2,Event)~core_hbv_Ab,KM_data_primary_event)
formula_retroviral <- survfit(Surv(time_2 ,Event)~hiv_treatment_at_hcv_diag.,data=KM_data_primary_event)
formula_peak_bil <- survfit(Surv(time_2,Event)~peak_bilirubin_binary, data=KM_data_primary_event,)
formula_alt<- survfit(Surv(time_2,Event)~Binary_peakALT,data=KM_data_primary_event)
formula_cd4 <- survfit(Surv(time_2,Event)~Group_CD4_count,data=KM_data_primary_event)
formula_cocaine <- survfit(Surv(time_2,Event)~Cocaine_use,data=KM_data_primary_event)
formula_heroin <- survfit(Surv(time_2,Event)~Heroin_use,data=KM_data_primary_event)
formula_meth <- survfit(Surv(time_2,Event)~Methamphetamine_use,data=KM_data_primary_event)
formula_gt <- survfit(Surv(time_2,Event)~combined_genotype,data=KM_data_primary_event)
formula_il28b <- survfit(Surv(time_2,Event)~IL28B_genotype_binary,data=KM_data_primary_event)
fits <- list(Age=formula_age,Gender=formula_gender,HIV=formula_hiv,
             PWID=formula_pwid,`Hepatitis B virus`=formula_hbv,`antiretroviral therapy`=formula_retroviral,
             `Peak Bilirubin`=formula_peak_bil,`Peak ALT`=formula_alt,`CD4 count`=formula_cd4,
             Cocaine=formula_cocaine,Heroin=formula_heroin,Methamphetamine=formula_meth,
             `Genotype`=formula_gt,IL28B=formula_il28b)

legend_titles <- list("Age","Gender","HIV status","PWID","Core HBV Ab",
                      "Patients on ARVS","peak Bilirubin","peak ALT","CD4 count",
                      "Cocaine Use","Heroin Use","Methamphetamine Use","Genotype",
                      "IL28B")
legend_labs <- list(c("\u2264 45",">45"),c("Male","Female"),c("hiv positive","hiv negative"),c("No","Yes"),c("Negative","Positive"),
                    c("No","Yes"),c("<20","\u2265 20"),c("\u2264 1000",">1000"),
                    c("\u2265 200","<200"),c("No","Yes"),
                    c("No","Yes"),c("No","Yes"),c("gt1","gt3","gt4","multiple genotypes"),c("CT/TT","CC"))

list_plots <- ggsurvplot_list(fit=fits,data=KM_data_primary_event,legend.title = legend_titles,legend.labs = legend_labs,
                              conf.int=TRUE,pval=TRUE,ylab="Proportion not clearing HCV",xlabs="Days")
#Indivudal ones
lapply(names(list_plots),
       function(x)ggsave(path="Output/survival analysis_uni_sc/",filename=paste(x,".pdf",sep=""),height = 5.83,width=8.27, units= "in",
                         plot=print(list_plots[[x]]),
                         device = cairo_pdf))
#Only plots
plots_only <- lapply(list_plots,function(x) x$plot)

#Arranged in a 2x2 plots but had to rendered manuually as symbols were not coming up properly
ggarrange(plotlist = plots_only[1:4],nrow=2,ncol=2,labels=c("A","B","C","D"))
ggarrange(plotlist = plots_only[5:8],nrow = 2,ncol=2,labels = c("E","F","G","H"))
ggarrange(plotlist = plots_only[9:12],nrow = 2,ncol=2,labels = c("I","J","K","L"))
ggarrange(plotlist = plots_only[13:14],nrow = 2,ncol=2,labels = c("M","N"))


#Only the significant KM ones####
significant_ones <- arrange_ggsurvplots(list_plots[c(3,7,11,13)],print = FALSE,ncol = 2,nrow=2)
ggsave("Output/survival analysis_uni_sc/significant_ones_KM.pdf",significant_ones,device = cairo_pdf,width = 297,height=210,units = "mm")
    


