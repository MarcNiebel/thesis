#Created by Marc Niebel April 2021
#Purpose is to carry out cox proportional hazards model on the acute cohort
#at a multivariable level

#Libraries needed
library(dplyr)
library(survival)
library(survminer)
library(gridExtra)

#Variables and time_1,time2 to event
variable_data <- read.csv("Output/survival analysis_uni_sc/survival_df.csv")

#Cleaning up what  is not needed
variable_data_cleaned <- variable_data  %>%
    select(-c(reinfection,spont_clearance,multiple_genotypes,acute_infection))

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

#Change name of column
variable_data_cleaned <- rename(variable_data_cleaned,HIV=hiv_pos)
variable_data_cleaned <-rename(variable_data_cleaned,`Peak Bilirubin`=peak_bilirubin_binary)
variable_data_cleaned<-rename(variable_data_cleaned,`Peak ALT`=Binary_peakALT)
variable_data_cleaned <- rename(variable_data_cleaned,Genotype=combined_genotype)
variable_data_cleaned <- rename(variable_data_cleaned,`Core HBV Antibody`=core_hbv_Ab)
variable_data_cleaned <- rename(variable_data_cleaned,`Heroin Use`=Heroin_use)
variable_data_cleaned<- rename(variable_data_cleaned,`Cocaine Use`= Cocaine_use)
variable_data_cleaned <- rename(variable_data_cleaned,`Methamphetamine Use`=Methamphetamine_use)

#Combining IL28B levels
variable_data_cleaned <- variable_data_cleaned %>% 
    mutate(IL28B_genotype_binary=case_when(IL28B_genotype == "CT"| 
                                               IL28B_genotype == "TT" ~"non-CC",
                                           IL28B_genotype =="CC" ~ "CC",
                                           TRUE ~ NA_character_)) %>%
    select(-IL28B_genotype)
variable_data_cleaned$IL28B_genotype_binary <- relevel(factor(variable_data_cleaned$IL28B_genotype_binary),ref = "non-CC")
variable_data_cleaned <- rename(variable_data_cleaned,IL28B=IL28B_genotype_binary)

#HIV positive only
variable_data_cleaned$Group_CD4_count <-relevel(factor(variable_data_cleaned$Group_CD4_count),ref=">=200")
variable_data_cleaned$hiv_treatment_at_hcv_diag. <- relevel(variable_data_cleaned$hiv_treatment_at_hcv_diag.,ref="No")
variable_data_cleaned <- rename(variable_data_cleaned, `CD4 count`=Group_CD4_count)
variable_data_cleaned <- rename(variable_data_cleaned,`Antiretroviral therapy`=hiv_treatment_at_hcv_diag.)

#Primary analysis####
#Model 1(AG)
model_1 <- coxph(Surv(time_1,time_2,Event)~Age+Gender+HIV+`Peak Bilirubin`+
                     `Peak ALT`+PWID+Genotype+cluster(record_id),data=variable_data_cleaned)
model_1_summary <- summary(model_1)
sink("Output/Multivariable_sc/primary_analysis_model_1.txt")
print(model_1_summary)
sink(file=NULL)

pdf(file="Output/Multivariable_sc/Forest_plot_primary_analysis_model_1.pdf", onefile = FALSE)
#Forest plot to summarise data
ggforest(model_1,data=variable_data_cleaned,main="Model 1")
dev.off()

#Exploratory analysis
#Model 2(AG)
#Added variables are recreational drug usage, Core HBV antibody and IL28B
#Removed is infected genotype as patients who spontaneously clear 
#cannot be typed and differential drug usage
model_2 <- update(model_1,~.+`Core HBV Antibody`+IL28B+`Cocaine Use`+
                      `Methamphetamine Use`+`Heroin Use`
                  -Genotype -PWID)
model_2_summary <- summary(model_2)
sink("Output/Multivariable_sc/exploratory_analysis_model_2.txt")
print(model_2_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/Forest_plot_exploratory_analysis_model_2.pdf",onefile = FALSE)
ggforest(model_2,data=variable_data_cleaned,main="Model 2")
dev.off()

#HIV positive patients####(AG)
hiv_positive_patients <- variable_data_cleaned %>% filter(HIV=="Yes")
model_3 <- coxph(Surv(time_1,time_2,Event)~Age+Gender+`Peak Bilirubin`+
                     `Peak ALT`+PWID+Genotype+`CD4 count`+
                     `Antiretroviral therapy`+cluster(record_id),data=hiv_positive_patients)
model_3_summary <- summary(model_3)
sink("Output/Multivariable_sc/hiv_positive_patients_only_model_3.txt")
print(model_3_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/Forest_plot_hiv_positive_patients_only_model_3.pdf",onefile = FALSE)
ggforest(model_3,data=hiv_positive_patients,main="Model 3")
dev.off()

#PWP model(total time)(model_1)
variable_data_cleaned_enum <- variable_data_cleaned %>%
    group_by(record_id) %>%
    mutate(enum=row_number())
model_4 <- coxph(Surv(time_1,time_2,Event)~Age+Gender+HIV+`Peak Bilirubin`+
                     `Peak ALT`+PWID+Genotype+cluster(record_id)+strata(enum),data=variable_data_cleaned_enum)
model_4_summary <-summary(model_4)
sink("Output/Multivariable_sc/pwp_total_time_model_4.txt")
print(model_4_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/pwp_total_time_model_4.pdf",onefile = FALSE)
ggforest(model_4,data=variable_data_cleaned_enum,main="Model 4")
dev.off()

#PWP model(gap-time)(model_1)
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

model_5 <- coxph(Surv(time_1,time_2,Event)~Age+Gender+HIV+`Peak Bilirubin`+
                     `Peak ALT`+PWID+Genotype+cluster(record_id)+strata(enum),data=variable_data_cleaned_diff_enum)
model_5_summary <- summary(model_5)
sink("Output/Multivariable_sc/pwp_gap_time_model_5.txt")
print(model_5_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/pwp_gap_time_model_5.pdf",onefile = FALSE)
ggforest(model_5,data=variable_data_cleaned_diff_enum,main="Model 1")
dev.off()

#Model 5 extended for exploratory(aka model_2)
model_6 <- update(model_5,~.+`Core HBV Antibody`+IL28B+`Cocaine Use`+
                      `Methamphetamine Use`+`Heroin Use`
                  -Genotype -PWID)
model_6_summary <- summary(model_6)
sink("Output/Multivariable_sc/pwp_gap_time_model_6.txt")
print(model_6_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/pwp_gap_time_model_6.pdf",onefile = FALSE)
ggforest(model_6,data=variable_data_cleaned_diff_enum,main="Model 2")
dev.off()

#Model 5  extended for HIV positive patients only(aka model_3)
hiv_positive_patients_enum <- variable_data_cleaned_diff_enum %>% filter(HIV=="Yes")
model_7 <- coxph(Surv(time_1,time_2,Event)~Age+Gender+`Peak Bilirubin`+
                     `Peak ALT`+PWID+Genotype+`CD4 count`+
                     `Antiretroviral therapy`+cluster(record_id)+strata(enum),data=hiv_positive_patients_enum)
model_7_summary <- summary(model_7)
sink("Output/Multivariable_sc/pwp_gap_time_model_7_hiv.txt")
print(model_7_summary)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/pwp_gap_time_model_7_hiv.pdf",onefile = FALSE)
ggforest(model_7,data=hiv_positive_patients_enum,main="Model 3")
dev.off()

#Proportional hazards assumption test####
#Using Schoenfeld residuals which should be independent of time
sch_residuals_model_1 <- cox.zph(model_1)
sch_residuals_model_2 <- cox.zph(model_2)
sch_residuals_model_3 <- cox.zph(model_3)
sch_residuals_model_4 <- cox.zph(model_4)
sch_residuals_model_5 <- cox.zph(model_5)
sch_residuals_model_6 <- cox.zph(model_6)
sch_residuals_model_7 <- cox.zph(model_7)

#Graphical diagnostics using ggcoxzph function which includes schoenfeld results
graphical_model_1 <- ggcoxzph(sch_residuals_model_1)
graphical_model_2 <- ggcoxzph(sch_residuals_model_2)
graphical_model_3 <- ggcoxzph(sch_residuals_model_3)
graphical_model_4 <- ggcoxzph(sch_residuals_model_4)
graphical_model_5 <- ggcoxzph(sch_residuals_model_5)
graphical_model_6 <- ggcoxzph(sch_residuals_model_6)
graphical_model_7 <- ggcoxzph(sch_residuals_model_7)

#Outputs of all seven
model_1_variable_on_different_pages <- marrangeGrob(graphical_model_1,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_1_sch_residuals.pdf",model_1_variable_on_different_pages)
model_2_variable_on_different_pages <- marrangeGrob(graphical_model_2,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_2_sch_residuals.pdf",model_2_variable_on_different_pages)
model_3_variable_on_different_pages <- marrangeGrob(graphical_model_3,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_3_sch_residuals.pdf",model_3_variable_on_different_pages)
model_4_variable_on_different_pages <- marrangeGrob(graphical_model_4,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_4_sch_residuals.pdf",model_4_variable_on_different_pages)
model_5_variable_on_different_pages <- marrangeGrob(graphical_model_5,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_5_sch_residuals.pdf",model_5_variable_on_different_pages)
model_6_variable_on_different_pages <- marrangeGrob(graphical_model_6,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_6_sch_residuals.pdf",model_6_variable_on_different_pages)
model_7_variable_on_different_pages <- marrangeGrob(graphical_model_7,ncol = 1,nrow = 1)
ggsave("Output/Multivariable_sc/model_7_sch_residuals.pdf",model_7_variable_on_different_pages)

#Changes need to be made to the ggforest function####
#https://stackoverflow.com/questions/63821274/plotting-a-cox-ph-model-using-ggforest-in-rstudio-when-a-factor-is-stratified
#Step 1: 
trace(ggforest, edit = TRUE)
#Replace at line 10-25 the function below and save. Run ggforest with model then.
allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if(var %in% colnames(data)) {
        if (terms[i] %in% c("factor", "character")) {
            adf <- as.data.frame(table(data[, var]))
            cbind(var = var, adf, pos = 1:nrow(adf))
        }
        else if (terms[i] == "numeric") {
            data.frame(var = var, Var1 = "", Freq = nrow(data), 
                       pos = 1)
        }
        else {
            vars = grep(paste0("^", var, "*."), coef$term, 
                        value = TRUE)
            data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                       pos = seq_along(vars))
        }
    } else {
        message(var, "is not found in data columns, and will be skipped.")
    }
    
})