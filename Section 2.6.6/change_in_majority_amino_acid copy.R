#Created by Marc Niebel Dec 2021
#Purpose was to look at the amino acid changes occuring from baseline

#Libraries needed
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)

#Path to excel data
excel_data <- "Data/Longitudinal.xlsx"
#HCV gene demarcation
gene_demarcation <- read_csv("Data/HCV_gene_demarcation.csv",col_types = cols(.default = "c"))
gene_demarcation <- gene_demarcation %>%
    mutate(Amino_acid=trimws(Amino_acid,which = "left"))
tab_names <- excel_sheets(path = excel_data)
#Putting all excel sheets together
list_all <- lapply(setNames(tab_names,tab_names),function(x)  read_excel(path=excel_data,sheet = x,col_types = "text"))

#Filter out ones which have no amino acid changes for later
no_amino_acid_changes <- lapply(list_all,function(x) dplyr::filter(x,amino_acid_changes=="None found"))
no_amino_acid_changes <- Filter(function(x) dim(x)[1]>0,no_amino_acid_changes)

#Filter ones with data
amino_acid_changes <- lapply(list_all,function(x) dplyr::filter(x,amino_acid_changes!="None found"))
amino_acid_changes <- Filter(function(x) dim(x)[1]>0,amino_acid_changes)
list_all_add_presence_ab <- lapply(amino_acid_changes,transform,pres_ab="1")

#Removing white space
removing_all_white_space  <- list_all_add_presence_ab %>% map(~map_df(.,~trimws(.,which = "left")))
#Binding them together
df_of_aa_changes <- bind_rows(removing_all_white_space,.id = "Patient")
#Still a problem so extracting number
df_of_aa_changes$amino_acid_changes <- as.numeric(sub("\\D+","",df_of_aa_changes$amino_acid_changes))
df_of_aa_changes$amino_acid_changes  <- as.character(df_of_aa_changes$amino_acid_changes)

#split the dataframes
df_of_aa_changes_split <- split(df_of_aa_changes,df_of_aa_changes$Patient)

#Adding on genes onto each  dataframe
adding_gene_demarcation <- lapply(df_of_aa_changes_split,
       right_join,gene_demarcation,by=c("amino_acid_changes"="Amino_acid"))
#Filling the NA in patient for each patient
filling_na_in_patient <- lapply(adding_gene_demarcation,fill,Patient)

#need to sort amino acid changes by number
mutate_amino_acids_values <- map(filling_na_in_patient,.%>% mutate_at("amino_acid_changes",as.numeric))
sorting_by_amino_acid_numbering <- lapply(mutate_amino_acids_values,function(df){arrange(df,amino_acid_changes)})

#Patients with no amino acid changes
adding_on_gene_dem_to_no_change_patients <- lapply(no_amino_acid_changes,function(x) cbind(x,gene_demarcation))
combine_no_change_patients<-bind_rows(adding_on_gene_dem_to_no_change_patients,.id = "Patient")
change_entry_combine_no_change_patients <- combine_no_change_patients %>%
    mutate_at("amino_acid_changes",str_replace,"None found","0") %>%
    rename(.,pres_ab=amino_acid_changes) %>%
    select(Patient,Amino_acid,pres_ab,Gene)

#Combine all data
combining_all_patients <- bind_rows(sorting_by_amino_acid_numbering,.id = "Patient")
colnames(combining_all_patients)[2] <-"Amino_acid"
combining_all_patients$Amino_acid <-as.character(combining_all_patients$Amino_acid)

#Filling in zeros
combining_all_patients$pres_ab[is.na(combining_all_patients$pres_ab)]<-0

#Bind both no change and changes together
all_data_inclusive <- bind_rows(combining_all_patients,change_entry_combine_no_change_patients)

relevant_cols <- all_data_inclusive %>%
    select(-c(Gene))
relevant_cols_wider <- relevant_cols %>%
    pivot_wider(names_from = "Amino_acid",values_from = "pres_ab")
relevant_cols_wider[] <-lapply(relevant_cols_wider,as.factor)
    
#Function to carry out shannon entropy alignment
entropy_column <- function(col){
    frequencies <-table(col)/length(col)
    df<-as.data.frame(frequencies)[,2]
    df<-df[df>0]
    -sum(df * log2(df))
}
#Entropy on each column of the  dataframe    
list_of_entropies <- apply(relevant_cols_wider[,-c(1)],2,entropy_column)
#Making it into a dataframe
df_of_entropies <-as.data.frame(list_of_entropies)
df_of_entropies$Amino_acid <- row.names(df_of_entropies)
#Adding on genes
df_of_entropies_genome_demaracation <- left_join(df_of_entropies,gene_demarcation,by="Amino_acid")
#Ensuring windows and genes are in correct order for plot
df_of_entropies_genome_demaracation$Amino_acid <- factor(df_of_entropies_genome_demaracation$Amino_acid,levels=df_of_entropies_genome_demaracation$Amino_acid)
df_of_entropies_genome_demaracation$Gene <- factor(df_of_entropies_genome_demaracation$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

#Plot across the genome
across_the_genome_plot <- ggplot(df_of_entropies_genome_demaracation,aes(x=Amino_acid,y=list_of_entropies,color=Gene,group=1))+geom_line()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key=element_blank(),
          axis.line=element_line(color="black"))+
    scale_y_continuous(expand=c(0,0),limits = c(0,1.00))+
    scale_color_npg()+labs(x="Amino acids across the genome",y="Shannon entropy")

list_of_labels_wanted <- c("384","412","421","460","485","519","536",
                           "570","581","597","645","746")
#Only E2
e2_only <- df_of_entropies_genome_demaracation %>% 
    filter(Gene=="E2")
e2_plot <- ggplot(e2_only,aes(x=Amino_acid,y=list_of_entropies,color=Gene,group=1))+geom_point(col="#00A087FF")+
    theme(legend.key=element_blank())+
    scale_x_discrete(breaks=list_of_labels_wanted,expand =c(0.01,0))+
    labs(x="Amino acids across the E2 region",y="Shannon entropy")+theme_classic()

#Counting the number of amino acids where there is change
counts_of_types  <- all_data_inclusive  %>%
    mutate_at(vars(pres_ab),as.numeric)
#Going to remove anything below 10%
count_of_amino_acids  <- counts_of_types %>%
    group_by(Amino_acid,Gene) %>%
    summarise(prop=sum(pres_ab)/n_distinct(Patient)) %>%
    arrange(desc(prop)) %>%
    filter(prop>0.1)
#Make a factor for plotting
count_of_amino_acids$Amino_acid <- factor(count_of_amino_acids$Amino_acid,levels=unique(count_of_amino_acids$Amino_acid[order(count_of_amino_acids$prop)]))
count_of_amino_acids$Gene <- factor(count_of_amino_acids$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
colors_needed <- c("E1"="#4DBBD5FF","E2"="#00A087FF","P7"="#3C5488FF","NS2"="#F39B7FFF","NS3"="#8491B4FF","NS4A"="#91D1C2FF","NS4B"="#DC0000FF","NS5A"="#7E6148FF","NS5B"="#B09C85FF")
table_of_targets <- ggplot(data=count_of_amino_acids,aes(y=prop,x=Amino_acid,fill=Gene))+
    geom_col()+coord_flip()+scale_fill_manual(values = colors_needed)+theme_classic()+labs(y="Proportion of sequences with amino acid change",
                                                                                       x="Amino acid position")

all_plots_together <- (across_the_genome_plot+table_of_targets)/e2_plot+plot_layout(heights = c(2,1))+plot_annotation(tag_levels = 'A')  & theme(plot.tag = element_text(face = "bold"))

ggsave("Output/plots_for_amino_acid_changes.pdf",all_plots_together,width = 11.69,height = 8.27)

#Data needed for downstream analysis
write_csv(relevant_cols_wider,"Output/amino_acid_changes_compared_to_baseline.csv")
patients <- relevant_cols_wider %>%
    select(Patient) %>%
    separate(Patient,c("Patient_new")) %>%
    distinct()

hla_data <- read_csv("Data/Baseline diversity/filtered_HLA_patients_presence_absence.csv")

relevant_patient_hla_data <- left_join(patients,hla_data,by=c("Patient_new"="record_id"))

write.csv(relevant_patient_hla_data,"Output/amino_acid_change_hla_data.csv")
