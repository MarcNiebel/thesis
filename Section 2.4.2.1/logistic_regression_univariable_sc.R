#Created by Marc Niebel Feb 2021
#Purpose is to create table of logistical regression results 
#with 95% CI and p-values

#Libraries needed
library(dplyr)
library(ggplot2)
library(jtools)

#Reading in data####
data <- read.csv("Output/univariable_sc/transformed_acute_data.csv")

#Preparing data for logistical regression####
# Ensuring the dependent variable is interogated as 0/1
converting_dependent_variable <- data %>%
    mutate(spont_clearance=case_when(spont_clearance == "Yes" ~ 1,
                                     TRUE ~ 0),
           spont_clearance=as.integer(spont_clearance))

#Change of reference for logistical regression to make sense
converting_dependent_variable$Gender <- relevel(converting_dependent_variable$Gender,ref = "Male")
converting_dependent_variable$hiv_pos <- relevel(converting_dependent_variable$hiv_pos,ref = "Yes")
converting_dependent_variable$core_hbv_Ab <-relevel(converting_dependent_variable$core_hbv_Ab,ref="Negative")
converting_dependent_variable$Cocaine_use <- relevel(converting_dependent_variable$Cocaine_use,ref = "No")
converting_dependent_variable$Heroin_use <- relevel(converting_dependent_variable$Heroin_use,ref="No")
converting_dependent_variable$Methamphetamine_use <- relevel(converting_dependent_variable$Methamphetamine_use,ref="No")
converting_dependent_variable$PWID <- relevel(converting_dependent_variable$PWID,ref = "No")
converting_dependent_variable$MSM <- relevel(converting_dependent_variable$MSM,ref = "Yes")

#HIV positive only
converting_dependent_variable$Group_CD4_count <-relevel(factor(converting_dependent_variable$Group_CD4_count),ref=">=200")
converting_dependent_variable$hiv_treatment_at_hcv_diag. <- relevel(converting_dependent_variable$hiv_treatment_at_hcv_diag.,ref="No")

#IL28B transformation
#Combine CT and TT due to low numbers of TT and not 
#wanting to use CC(associated with SC) as reference
combined_IL28B_levels <- converting_dependent_variable %>% 
    mutate(IL28B_genotype_binary=case_when(IL28B_genotype == "CT"| 
                                           IL28B_genotype == "TT" ~"non-CC",
                                           IL28B_genotype =="CC" ~ "CC",
                                     TRUE ~ NA_character_)) %>%
    select(-IL28B_genotype)

IL28B_new <- combined_IL28B_levels %>%
    select(IL28B_genotype_binary,spont_clearance) %>%
    mutate(spont_clearance=factor(spont_clearance,levels=c(1,0)))
xlabels <- c("Yes","No")
pdf("Output/univariable_sc/IL28B_binary.pdf")
ggplot(IL28B_new,aes(x=spont_clearance,fill=IL28B_genotype_binary))+
    xlab("Spontaneous Clearance")+
    ylim(0,205)+
    geom_bar(position = "dodge")+
    geom_text(aes(label = ..count..), 
              stat = "count", vjust = -0.5, 
              colour = "black",
              position = position_dodge(width = 1)) +
    labs(fill="IL28B designation")+
    theme_classic()+
    scale_x_discrete(labels=xlabels)
dev.off()
    
#Making non-CC reference
combined_IL28B_levels$IL28B_genotype_binary <- relevel(factor(combined_IL28B_levels$IL28B_genotype_binary),ref = "non-CC")

#Dataframe which can be used in downstream analysis
write.csv(combined_IL28B_levels,"Output/logistical_regression_univariable_sc/all_data_transformed.csv",row.names = FALSE)

#function to carry out general linearised model for each variable
logistical_regression <- function(dataframe){
    #Empty dataframe
    results_log_regression <- data.frame()
    #Variables of interest
    variables_of_interest <- dataframe[-c(1:3)]
    for(i in seq_along(variables_of_interest)){
        #Generate formula for each variable
        gen_formula <- as.formula(paste("spont_clearance ~",variables_of_interest[i]))
        #GLM model
        log_reg <- glm(gen_formula,data=dataframe,family="binomial")
        #Summary statistic and name
        variable=names(variables_of_interest[i])
        summary_log_reg <-summary(log_reg)
        odds_ratio <- exp(coef(summary_log_reg))[2]
        p_value <- coef(summary_log_reg)[8]
        confint_lower <- exp(confint(log_reg))[2]
        confint_higher <- exp(confint(log_reg))[4]
        #Temporary dataframe
        temp_df <- data.frame(Variables=variable,
                              OddsRatio=odds_ratio,p.value = p_value,
                              CI_lower=confint_lower,CI_higher=confint_higher)
        #Binding of all the data together
        results_log_regression <- rbind(results_log_regression,temp_df)
    }
    return(results_log_regression)
}

#Note the above function can only be used for an independent variable with
# two levels e.g. Yes/No or continous variable
#Therefore infected genotype will be accessed separately
infected_genotype_separated <- combined_IL28B_levels[-c(18)]
#Logistical regression of bivariate levels and age
two_level_logistical_regression <- logistical_regression(infected_genotype_separated)

#Rounding the numbers to 4 decimal places
two_level_logistical_regression[,-1] <- round(two_level_logistical_regression[,-1],4)

write.csv(two_level_logistical_regression,file="Output/logistical_regression_univariable_sc/logistical_regression_summary.csv")

#Multiple level variables
infected_genotype_variables <- combined_IL28B_levels[c(1:3,18)]

#Infected genotype
infected_genotype_fit <-glm(spont_clearance ~ Genotypes, 
                        data=infected_genotype_variables,family="binomial")
logit_genotype_co <-summ(infected_genotype_fit,exp=TRUE,digits = 4)
sink("Output/logistical_regression_univariable_sc/Infected_genotype_regression_summary.txt")
print(logit_genotype_co)
sink(file=NULL)
