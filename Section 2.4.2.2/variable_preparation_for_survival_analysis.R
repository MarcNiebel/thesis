#Created by Marc Niebel March 2021
#Purpose is to prepare data for survival analysis

#libraries needed
library(dplyr)
library(tidyr)

#Preparing data####
#Sourcing the dataframe from the data folder
source("./Data/AcuteHCV_R_2021-04-22_1259.r", chdir = TRUE)

#Added variables(peak bilirubin and peak ALT) of interest are for re-infection patients
#Have added more than one genotype and second clinical genotype as well as not done in first analysis
variables_of_interest <- data %>% 
    select(record_id,acute.factor,sc.factor,gender.factor,age,hiv.factor,bil_peak_1,
           alt_peak,risk___1.factor,hbvsag_pcr.factor,cab.factor,hiv_tx_at_hcv_diagnosis.factor,
           cd4_at_hcv_diagnosis,ifnl4_860.factor,risk___11.factor,risk___10.factor,risk___2.factor,
           drugs___10.factor:drugs___13.factor,drugs___1.factor:drugs___3.factor,drugs___20.factor,
           clinical_genotype___1.factor:clinical_genotype___4.factor,reinfection.factor,alt_peak_2,
           bil_peak_2,more_than_one_genotype.factor,clinical_genotype_2___1.factor,clinical_genotype_2___2.factor,
           clinical_genotype_2___3.factor,clinical_genotype_2___4.factor)

#Removed chronic patients####
removed_chronic_patients <- variables_of_interest %>% 
    filter(acute.factor == "Yes")

#Change names of columns####
colnames(removed_chronic_patients) <- c("record_id","acute_infection","spont_clearance",
                                        "Gender","Age","hiv_pos","peak_bilirubin",
                                        "peak_ALT","PWID","chronic_HBV","core_hbv_Ab",
                                        "hiv_treatment_at_hcv_diag.","CD4_count_at_hcv_diag.",
                                        "IL28B_genotype","MSM_unknown","MSM_insertive",
                                        "MSM_receptive","Cocaine_smoked","Cocaine_injected",
                                        "Cocaine_intranasal","Cocaine_rectal","Methamphetamine_smoked",
                                        "Methamphetamine_iv","Methamphetamine_other","Heroin_iv",
                                        "clinical_gt1","clinical_gt2","clinical_gt3","clinical_gt4",
                                        "reinfection","peak_ALT_2","peak_bilirubin_2","multiple_genotypes",
                                        "clinical_gt1_2","clinical_gt2_2","clinical_gt3_2","clinical_gt4_2")

#Transformations of data####
#Peak_bilirubin
#Using clinical cut-off of 20 mg/dL
assign_peak_bil_group <-removed_chronic_patients %>% 
    mutate(peak_bilirubin_binary=case_when(peak_bilirubin >=20 ~">= 20",
                                           peak_bilirubin < 20 ~"< 20",
                                           TRUE ~ NA_character_),
           peak_bilirubin_binary=as.factor(peak_bilirubin_binary)) %>%
    select(-peak_bilirubin)

assign_peak_bil_2_group <-assign_peak_bil_group %>% 
    mutate(peak_bilirubin_binary_2=case_when(peak_bilirubin_2 >=20 ~">= 20",
                                           peak_bilirubin_2 < 20 ~"< 20",
                                           TRUE ~ NA_character_),
           peak_bilirubin_binary_2=as.factor(peak_bilirubin_binary_2)) %>%
    select(-peak_bilirubin_2)

#Peak ALT
#Using clinical cut-off of 1000 IU/ml
assign_peak_alt_group <- assign_peak_bil_2_group %>%
    mutate(Binary_peakALT=case_when(peak_ALT >=1000 ~ ">=1000",
                                    peak_ALT < 1000 ~"<1000",
                                    TRUE ~ NA_character_),
           Binary_peakALT =as.factor(Binary_peakALT)) %>%
    select(-peak_ALT)

assign_peak_alt_2_group <- assign_peak_alt_group %>%
    mutate(Binary_peakALT_2=case_when(peak_ALT_2 >=1000 ~ ">=1000",
                                    peak_ALT_2 < 1000 ~"<1000",
                                    TRUE ~ NA_character_),
           Binary_peakALT_2 =as.factor(Binary_peakALT_2)) %>%
    select(-peak_ALT_2)

#CD4 count
#Using clinical cut-off of 200 cells/mm^3(severally immunocompromised)
assign_cd4_group <- assign_peak_alt_2_group %>%
    mutate(Group_CD4_count=case_when(CD4_count_at_hcv_diag. < 200 ~"<200",
                                     CD4_count_at_hcv_diag. >=200 ~">=200",
                                     TRUE ~ NA_character_),
           Group_CD4_count=as.factor(Group_CD4_count)) %>%
    select(-CD4_count_at_hcv_diag.)

#MSM
#Due to low numbers of the different transmission routes MSM will be amalgamated
MSM_risk <- assign_cd4_group %>% 
    mutate(MSM=case_when(MSM_unknown =="Checked"|
                             MSM_insertive =="Checked"|
                             MSM_receptive=="Checked" ~ "Yes",
                         TRUE ~ "No"),
           MSM=as.factor(MSM)) %>%
    select(-MSM_unknown,-MSM_receptive,-MSM_insertive)

#Cocaine Usage
#Low numbers in each therefore amalgamate
Cocaine_usage <- MSM_risk %>% 
    mutate(Cocaine_use=case_when(Cocaine_smoked=="Checked"|
                                     Cocaine_injected=="Checked"|
                                     Cocaine_intranasal =="Checked"|
                                     Cocaine_smoked == "Checked" ~ "Yes",
                                 TRUE ~"No"),
           Cocaine_use=as.factor(Cocaine_use))%>%
    select(-Cocaine_injected,-Cocaine_smoked,-Cocaine_rectal,-Cocaine_intranasal)

#Heroin Usage
#Changing designation from heroin iv to heroin usage
Heroin_usage <- Cocaine_usage %>% 
    mutate(Heroin_use=case_when(Heroin_iv=="Checked" ~ "Yes",
                                TRUE ~"No"),
           Heroin_use=as.factor(Heroin_use)) %>%
    select(-Heroin_iv)

#Methamphetamine Usage
#Low numbers lead to amalgamation of data
Meth_usage <- Heroin_usage %>% 
    mutate(Methamphetamine_use=case_when(Methamphetamine_smoked =="Checked"|
                                             Methamphetamine_iv =="Checked"|
                                             Methamphetamine_other =="Checked" ~"Yes",
                                         TRUE ~ "No"),
           Methamphetamine_use=as.factor(Methamphetamine_use)) %>%
    select(-Methamphetamine_iv,-Methamphetamine_smoked,-Methamphetamine_other)

#Chronic HBV
#Chronic HBV positive has no spontaneous clearers so it will be removed
removal_chronic_hbv_variable <- Meth_usage %>% select(-chronic_HBV)

#PWID
#Changing checked/unchecked to Yes/No
removal_chronic_hbv_variable <-removal_chronic_hbv_variable %>% 
    mutate(PWID=case_when(PWID == "Checked" ~ "Yes",
                          TRUE ~ "No"))
#Clinical genotype
clinical_genotype <- removal_chronic_hbv_variable %>% 
    select(record_id,clinical_gt1:clinical_gt4)

#Change into long format
clinical_genotype_long_format <- clinical_genotype %>% 
    pivot_longer(cols = clinical_gt1:clinical_gt4,names_to = "clinical_genotype")

#Only one gt 2 so remove but need for summary table
clinical_genotype_long_format %>% 
    filter(value =="Checked") %>% 
    group_by(clinical_genotype) %>% 
    summarise(n=n())

#Only the infected genotype
clinical_genotype_only <- clinical_genotype_long_format %>% 
    filter(value =="Checked") %>%
    select(-value)
#Join the infected data to main dataframe
transformed_data <- left_join(removal_chronic_hbv_variable,
                              clinical_genotype_only,by="record_id") %>%
    select(-c(clinical_gt1:clinical_gt4))

#secondary genotype
secondary_genotype <- removal_chronic_hbv_variable %>% 
    select(record_id,clinical_gt1_2:clinical_gt4_2)

#Change into long format
secondary_genotype_long_format <- secondary_genotype %>% 
    pivot_longer(cols = clinical_gt1_2:clinical_gt4_2,names_to = "secondary_genotype")

#Summary table
secondary_genotype_long_format %>% 
    filter(value =="Checked") %>% 
    group_by(secondary_genotype) %>% 
    summarise(n=n())

#Only the secondary genotype
secondary_genotype_only <- secondary_genotype_long_format %>% 
    filter(value =="Checked") %>%
    select(-value)
#Join the secondary data to main dataframe
transformed_data_sec_infection_added <- left_join(transformed_data,
                              secondary_genotype_only,by="record_id") %>%
    select(-c(clinical_gt1_2:clinical_gt4_2))

genotype_data <- transformed_data_sec_infection_added %>%
    select(record_id,clinical_genotype,secondary_genotype)

#Modifying clinical genotype
modifing_column_genotype_data <-genotype_data %>% 
    separate(clinical_genotype,c(NA,"clinical_genotype"))
modifying_column_secondary_genotype_data <-modifing_column_genotype_data %>%
    separate(secondary_genotype,c(NA,"secondary_genotype",NA))

#Combining both genotype results together
combined_genotype <-modifying_column_secondary_genotype_data %>% 
    unite("combined_genotype", clinical_genotype:secondary_genotype, 
          na.rm = TRUE, remove = FALSE) %>%
    select(-c(clinical_genotype,secondary_genotype))

#Replace blank cells with NA
combined_genotype$combined_genotype[combined_genotype$combined_genotype==''] <- NA

#Small numbers for the multiple genotypes therefore combined
combined_genotype %>%
    group_by(combined_genotype) %>%
    summarise(n=n())

#Different genotypes combined. Note also gt1_gt1 needs to be mutated to gt1 only
different_genotypes_combined <- combined_genotype %>%
    mutate(combined_genotype=case_when(combined_genotype =="gt1_gt3"|
                          combined_genotype == "gt1_gt4"|
                          combined_genotype == "gt3_gt1"|
                          combined_genotype == "gt4_gt1" ~ "multiple_genotypes",
                          combined_genotype =="gt1_gt1" ~ "gt1",
                          TRUE ~ combined_genotype))
#Numbers in each group after manipulation
different_genotypes_combined %>%
    group_by(combined_genotype) %>%
    summarise(n=n())

#Adding that onto entire dataframe
transformed_data_combined_data_added <- 
    left_join(transformed_data_sec_infection_added,
              different_genotypes_combined,by="record_id") %>%
    select(-c(clinical_genotype,secondary_genotype))

#Re_infection patients
re_infection_patients <- transformed_data_combined_data_added %>%
    filter(reinfection=="Yes")
#making the data longer for peak bilirubin
re_infection_bil_longer_format <- pivot_longer(re_infection_patients,
                                               cols=peak_bilirubin_binary:peak_bilirubin_binary_2,
                                               names_to = "peak_bilirubin_binary")  %>%
    select(-c(peak_bilirubin_binary,Binary_peakALT,Binary_peakALT_2,Binary_peakALT,Binary_peakALT_2))
names(re_infection_bil_longer_format)[19] <- "peak_bilirubin_binary"

#making the data longer for peak ALT
re_infection_alt_longer_format <- pivot_longer(re_infection_patients,
                                               cols=Binary_peakALT:Binary_peakALT_2,
                                               names_to = "peak_alt_binary") %>% 
    select(record_id,value)
names(re_infection_alt_longer_format)[2] <- "Binary_peakALT"

#Peak ALT and bilirubin manipulation complete
re_infection_patients_complete <-bind_cols(re_infection_alt_longer_format,re_infection_bil_longer_format) %>%
    select(-record_id1)

#Removing re-infection patients from original dataset
transformed_data_combined_data_added$reinfection[is.na(transformed_data_combined_data_added$reinfection)] <- "No"
re_infection_patients_removed <- transformed_data_combined_data_added %>%
    filter(reinfection =="No") %>% 
    select(-c(peak_bilirubin_binary_2,Binary_peakALT_2))

#Read in the time to event .csv with the status of event
time_to_event <-read.csv("data/time_to_event_multiple_instances.csv")

#Adding on time to event for mono infected patients
single_infection_patients <- left_join(re_infection_patients_removed,time_to_event,by="record_id")
single_infection_patients <- single_infection_patients %>% select(-c(Reinfection))

#Pulling out re_infected patients
time_to_event_re_infected_patients <- time_to_event %>%
    filter(Reinfection =="Yes")

#Problem samples for re_infection data
problem_samples <- c("56","60","68","77","106","220","G50")

#Remove from re_infections_complete
problem_samples_removed_from_metadata <- re_infection_patients_complete %>% 
    filter (!record_id %in% problem_samples)

#Remove from time_to_event_re_infected_patients
problem_samples_removed_from_time_to_event <- time_to_event_re_infected_patients %>% 
    filter (!record_id %in% problem_samples)  %>% select(-Reinfection)

no_problem_samples_combined <- bind_cols(problem_samples_removed_from_metadata,
                                         problem_samples_removed_from_time_to_event) %>%
    select(-record_id1)

problem_samples_from_metadata <- re_infection_patients_complete %>% 
    filter (record_id %in% problem_samples)

#This is to do with P106 for which we only have two measurements so last one will be deleted
problem_samples_from_time_to_event <- time_to_event_re_infected_patients %>% 
    filter (record_id %in% problem_samples)
removed_last_event_of_patient <- problem_samples_from_time_to_event[-3,]
metadata_for_patient <- problem_samples_from_metadata %>% filter(record_id=="106")
patient_metadata_time_combined <-bind_cols(removed_last_event_of_patient,metadata_for_patient)

#Resolving remaining discrepancies which 
#pertain to not having the intervals for survival analysis
resolved_re_infection_patients <- problem_samples_from_metadata[-c(2,4,6,8,9,10,12,13),]
resolved_re_infection_patients_metadata <- left_join(resolved_re_infection_patients,time_to_event,by="record_id")

#All re_infections and discrepancies together
df_binding_re_infections_discrep <- bind_rows(list(no_problem_samples_combined,patient_metadata_time_combined,
                       resolved_re_infection_patients_metadata)) %>%  select(-c(Reinfection,record_id1))

all_patients_together <- bind_rows(single_infection_patients,df_binding_re_infections_discrep)

#Write out to csv
#Writing the data for downstream analysis####
#Replace genotype 2 result with NA
all_patients_together_removed_gt2 <- all_patients_together %>%
    mutate(combined_genotype=case_when(combined_genotype =="gt2" ~ NA_character_,
                                       TRUE ~ combined_genotype))
write.csv(all_patients_together_removed_gt2,file="Output/survival analysis_uni_sc/survival_df.csv",row.names = FALSE)

#Follow-up time summary
follow_up_time <- data %>% 
    filter(acute.factor=="Yes") %>%
    select(follow_up)

quantile(follow_up_time$follow_up)
    
