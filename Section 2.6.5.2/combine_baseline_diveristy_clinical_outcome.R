#Created by Marc Niebel Nov 2021
#The purpose of this script is to first combine all patient
#diversity windows, add clinical outcome and other metadata
#before carrying visualisations at both patient, gene and window levels
#and including statistics where possible

#Libraries needed
library(readr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(gtools)
library(imputeTS)
library(patchwork)
library(ggforce)
library(stringr)
library(tidyr)

#All the window data
patient_all_windows_files <-list.files(path = "Output/Baseline diversity",pattern="*window_diversity.csv",full.names = TRUE)

#The gini_simpson of all windows for all patients
gini_simpson_all_window_data <- lapply(patient_all_windows_files,read_csv)

#Bringing all files into one dataframe
all_patients_windows_together <- bind_rows(gini_simpson_all_window_data,.id = "type")

#Filtering out P13 secondary and P110 secondary which will be addressed later
P13_P110_secondary_filtered_out <- all_patients_windows_together  %>%
    filter(!(type %in% c("43","49")))

#Reading in the clinical outcome for the patients analysed
baseline_metadata <- read_csv("Data/Baseline diversity/metadata_baslines.csv")

#Remove P13 and P110 Secondary
baseline_metadata_P13_P110_removed <- baseline_metadata %>% filter(!(id  %in% c("P110_230210_S59_R1_001_190409",
                                          "P13_071204_input_S7S_P13_071204_input_S7R1_131201")))

all_patients_windows_together_co <-left_join(P13_P110_secondary_filtered_out,
                                             baseline_metadata_P13_P110_removed,
                                             by=c("patient_id"="patient"))
#P13 and P110 secondary 
P13_P110_secondary_filtered <- all_patients_windows_together  %>%
    filter(type %in% c("43","49"))

baseline_metadata_P13_P110 <- baseline_metadata %>% 
    filter(id  %in% c("P110_230210_S59_R1_001_190409",
                        "P13_071204_input_S7S_P13_071204_input_S7R1_131201"))
P13_P110_secondary_co <-left_join(P13_P110_secondary_filtered,
                                             baseline_metadata_P13_P110,
                                             by=c("patient_id"="patient"))
#Bind both datasets together
all_patients_together_including_secondary <- bind_rows(all_patients_windows_together_co,
                                                       P13_P110_secondary_co)
all_patients_together_including_secondary <- all_patients_together_including_secondary %>%
    filter(!(patient_id %in% c("P220","G39")))
#Primary infection only####
primary_infection_only <- all_patients_together_including_secondary %>%
    filter(type.y=="Primary")
#Mean of each patient
mean_of_each_patient <- primary_infection_only %>%
    group_by(patient_id,clinical_outcome) %>%
    summarise(mean_by_patient=mean(gini_simpson,na.rm = TRUE))

patient_boxplot <- ggboxplot(mean_of_each_patient, y="mean_by_patient",x="clinical_outcome",
          color="clinical_outcome",outlier.shape=NA,
          xlab = "Clinical Outcome",ylab="Average Gini-Simpson Index",legend="none",add = "jitter")+
    stat_compare_means(label="p.format",label.x = 1.5)

#Mean by region####
mean_by_region <- primary_infection_only %>%
    group_by(Gene,clinical_outcome,patient_id) %>%
    summarise(mean_by_region=mean(gini_simpson,na.rm=TRUE))

#Ensure genes are in correct order
mean_by_region$Gene <- factor(mean_by_region$Gene,
                                              levels = c("Core","E1","E2","P7","NS2",
                                                         "NS3","NS4A","NS4B","NS5A","NS5B"))
region_plot <- ggboxplot(mean_by_region, y="mean_by_region",x="Gene",
          color="clinical_outcome",add.params = list(alpha = 0.3),outlier.shape = NA,
          xlab = "Genome Regions",ylab="Average Gini-Simpson Index",legend.title="Clinical Outcome",add = "jitter")+
    stat_compare_means(aes(group=clinical_outcome),label="p.format")+
    ylim(0,1)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

region_plot_facet <- facet(region_plot,facet.by = "Gene",scales = "free_x",nrow = 2,ncol = 5)

#Bringing the two together i.e. patient and gene level
patient_gene_plot_together <- (patient_boxplot|region_plot_facet)+
    plot_annotation(tag_levels = 'A')  & theme(plot.tag = element_text(face = "bold"))
ggsave("Output/patient_gene_diversity.pdf",patient_gene_plot_together,width=11.69,height = 8.27)
#Mean of each window of all primary infected patients####
mean_of_each_window <- primary_infection_only %>%
    group_by(Window,Gene,clinical_outcome,.groups="keep") %>%
    summarise(mean_by_window=mean(gini_simpson,na.rm = TRUE))
#Ensuring the order of the windows is maintained
mean_of_each_window <- mean_of_each_window[mixedorder(mean_of_each_window$Window),]
mean_of_each_window$Window <- factor(mean_of_each_window$Window,
                                     levels=unique(mean_of_each_window$Window))

#Change the order of levels to ensure correct gene order
mean_of_each_window$Gene <- factor(mean_of_each_window$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

all_window_plot_together <- ggscatter(mean_of_each_window, "Window", "mean_by_window",
          shape = "clinical_outcome",color = "Gene",xlab="Windows across genome",
          ylab="Average Gini-Simpson Index")+scale_color_npg()+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank())+
    scale_x_discrete(expand=expansion(add=50))+
    guides(shape=FALSE)
across_the_genome_by_clinical_outcome_plot <- facet(all_window_plot_together,facet.by = "clinical_outcome",nrow = 2,ncol = 1)
ggsave("Output/window_summary_diversity.pdf",across_the_genome_by_clinical_outcome_plot,width=11.69,height = 8.27)

#Getting all diversity  plots across the genome for all patients####
diversity_data_primary <- primary_infection_only

#Ensuring the order of the windows is maintained
diversity_data_primary <- diversity_data_primary[mixedorder(diversity_data_primary$Window),]
diversity_data_primary$Window <- factor(diversity_data_primary$Window,
                                     levels=unique(diversity_data_primary$Window))

#Change the order of levels to ensure correct gene order
diversity_data_primary$Gene <- factor(diversity_data_primary$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

diversity_data_primary_split <- split(diversity_data_primary,diversity_data_primary$patient_id)

diversity_plots <-
     lapply(seq_along(diversity_data_primary_split),
            function(df){
                ggplot_na_distribution(diversity_data_primary_split[[df]][4],color_missing = "lightgrey",color_missing_border = "lightgrey",
                                       xlab="Windows across genome", ylab = "Gini-Simpson Index",color_lines = "white",
                                       title = NULL,subtitle = NULL,theme=theme_classic())+
                    ggplot2::geom_point(data=diversity_data_primary_split[[df]],ggplot2::aes(x=Window,y=gini_simpson,
                                                                                 color=Gene))+
                    ggplot2::ylim(0,1)+
                    ggplot2::theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.border = element_rect(color="black",fill = NA))+
                    scale_color_npg()+
                    scale_x_discrete(expand=expansion(add=50))+
                    facet_wrap(~patient_id+clinical_outcome)
            })
all_diversity_plots <-  ggarrange(plotlist = diversity_plots,nrow=4,ncol=3,common.legend = TRUE)

ggexport(all_diversity_plots,filename = "Output/all_diversity_plots_together.pdf",width = 8.27,height = 11.69)

G10_plot <- diversity_plots[[1]]
P4_plot <- diversity_plots[[81]]
plot_list<- list(G10_plot,P4_plot)
thesis_plots <- ggarrange(plotlist = plot_list,nrow = 1,ncol = 2,common.legend = TRUE,labels =c ("A","B"))
ggsave("Output/G10_P4_diveristy_across_genome.pdf",thesis_plots,width = 11.69,height = 8.27)
#Looking at the data longitudinally
#Patient
mean_of_each_patient_with_time <- primary_infection_only %>%
    group_by(patient_id,clinical_outcome,days_diff) %>%
    summarise(mean_by_patient=mean(gini_simpson,na.rm = TRUE))
scatterplot_over_time <- ggscatter(mean_of_each_patient_with_time,x="days_diff",
                                   y="mean_by_patient",ylab="mean Gini-Simpson Index",
                                   xlab = "Timespan between EDI and sample taken(days)",
                                   facet.by = "clinical_outcome",color = "clinical_outcome")+
    scale_color_discrete(name="Clinical Outcome")
#By gene
mean_of_each_gene_with_time <- primary_infection_only %>%
    group_by(patient_id,Gene,clinical_outcome,days_diff) %>%
    summarise(mean_by_gene=mean(gini_simpson,na.rm = TRUE))
mean_of_each_gene_with_time$Gene <- factor(mean_of_each_gene_with_time$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
scatterplot_over_time_gene <- ggscatter(mean_of_each_gene_with_time,x="days_diff",
                                        y="mean_by_gene",
                                        ylab="mean Gini-Simpson Index",
                                        xlab = "Timespan between EDI and sample taken(days)",
                                        color = "clinical_outcome")+
    scale_color_discrete(name="Clinical Outcome")
scatterplot_over_time_gene_facet <- facet(scatterplot_over_time_gene,facet.by = "Gene",ncol = 5,nrow = 2)

longitudinal_data_together <- (scatterplot_over_time|scatterplot_over_time_gene_facet)+plot_layout(widths = c(1,2),guides = 'collect')+plot_annotation(tag_levels = 'A')&
    theme(legend.position = 'top',plot.tag = element_text(face = "bold"))

ggsave("Output/continous_time_patient_gene.pdf",longitudinal_data_together,width = 11.69,height = 8.27)

#Secondary infection####
secondary_infection_only <- all_patients_together_including_secondary %>%
    filter(type.y=="Secondary")
#Looking at difference at patient level
#Mean of each patient
mean_of_each_patient_secondary <- secondary_infection_only %>%
    group_by(patient_id,clinical_outcome) %>%
    summarise(mean_by_patient=mean(gini_simpson,na.rm = TRUE))

patient_boxplot_secondary <- ggboxplot(mean_of_each_patient_secondary, y="mean_by_patient",x="clinical_outcome",
                             color="clinical_outcome",outlier.shape=NA,
                             xlab = "Clinical Outcome",ylab="Average Gini-Simpson Index",legend="none",add = "jitter")+
    stat_compare_means(label="p.format",label.x = 1.5)

#Mean by region####
mean_by_region_sec <- secondary_infection_only %>%
    group_by(Gene,clinical_outcome,patient_id) %>%
    summarise(mean_by_region=mean(gini_simpson,na.rm=TRUE))

#Ensure genes are in correct order
mean_by_region_sec$Gene <- factor(mean_by_region_sec$Gene,
                              levels = c("Core","E1","E2","P7","NS2",
                                         "NS3","NS4A","NS4B","NS5A","NS5B"))
region_plot_sec <- ggboxplot(mean_by_region_sec, y="mean_by_region",x="Gene",
                         color="clinical_outcome",add.params = list(alpha = 0.3),outlier.shape = NA,
                         xlab = "Genome Regions",ylab="Average Gini-Simpson Index",legend.title="Clinical Outcome",add = "jitter")+
    stat_compare_means(aes(group=clinical_outcome),label="p.format",label.y = 0.6)+
    ylim(0,0.75)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
region_plot_sec_facet <-facet(region_plot_sec,facet.by = "Gene",scales = "free_x",nrow = 2,ncol = 5)

#Looking at averages across all windows####
mean_of_each_window_secondary <- secondary_infection_only %>%
    group_by(Window,Gene,clinical_outcome,.groups="keep") %>%
    summarise(mean_by_window=mean(gini_simpson,na.rm = TRUE))
#Ensuring the order of the windows is maintained
mean_of_each_window_secondary <- mean_of_each_window_secondary[mixedorder(mean_of_each_window_secondary$Window),]
mean_of_each_window_secondary$Window <- factor(mean_of_each_window_secondary$Window,
                                               levels=unique(mean_of_each_window_secondary$Window))

#Change the order of levels to ensure correct gene order
mean_of_each_window_secondary$Gene <- factor(mean_of_each_window_secondary$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

all_window_plot_together_secondary <- ggscatter(mean_of_each_window_secondary, "Window", "mean_by_window",
                                                shape = "clinical_outcome",color = "Gene",xlab="Windows across genome",
                                                ylab="Average Gini-Simpson Index")+scale_color_npg()+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank())+
    scale_x_discrete(expand=expansion(add=50))+
    guides(shape=FALSE)
across_the_genome_by_clinical_outcome_plot_secondary <- facet(all_window_plot_together_secondary,facet.by = "clinical_outcome",nrow = 2,ncol = 1)

patch <- (patient_boxplot_secondary+region_plot_sec_facet)+plot_layout(widths = c(1,1.5))
secondary_plots_together <- patch/
    across_the_genome_by_clinical_outcome_plot_secondary+plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(face = "bold"))

ggsave("Output/secondary_plots_secondary.pdf",secondary_plots_together,width=8.27,height = 11.69)

#Look at primary versus reinfection
#Patient level
mean_of_each_patient$type <-"Primary"
mean_of_each_patient_secondary$type <- "Secondary"
both_primary_secondary_together <- bind_rows(mean_of_each_patient,mean_of_each_patient_secondary)
both_primary_secondary_together$type <- factor(both_primary_secondary_together$type)
#Again numbers are so small not sure I want
colour_palette_to_be_used <-c("#20854EFF","#7876B1FF")
primary_secondary_diversity <- ggboxplot(both_primary_secondary_together,y="mean_by_patient",x="type",
                                              add = "jitter", color ="type",
                                              add.params = list(alpha = 0.3),outlier.shape = NA,palette =colour_palette_to_be_used )+
    stat_compare_means(label="p.format",label.x = 1.5)+
    theme(legend.position = "none")+labs(y="Average Gini-Simpson Index")
primary_secondary_diversity_facet <- facet(primary_secondary_diversity,facet.by = "clinical_outcome")+theme(legend.position = "none")

#Both plots together
primary_sec_together <- primary_secondary_diversity/primary_secondary_diversity_facet+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold"))

ggsave("Output/primary_secondary_together.pdf",primary_sec_together,width = 8.27,height = 11.69)

#Look at diversity across the genome####
diversity_data_secondary <- secondary_infection_only

#Ensuring the order of the windows is maintained
diversity_data_secondary <- diversity_data_secondary[mixedorder(diversity_data_secondary$Window),]
diversity_data_secondary$Window <- factor(diversity_data_secondary$Window,
                                        levels=unique(diversity_data_secondary$Window))

#Change the order of levels to ensure correct gene order
diversity_data_secondary$Gene <- factor(diversity_data_secondary$Gene,levels=c("Core","E1","E2","P7","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

diversity_data_secondary_split <- split(diversity_data_secondary,diversity_data_secondary$patient_id)

diversity_plots_secondary <-
    lapply(seq_along(diversity_data_secondary_split),
           function(df){
               ggplot_na_distribution(diversity_data_secondary_split[[df]][4],color_missing = "lightgrey",color_missing_border = "lightgrey",
                                      xlab="Windows across genome", ylab = "Gini-Simpson Index",color_lines = "white",
                                      title = NULL,subtitle = NULL,theme=theme_classic())+
                   ggplot2::geom_point(data=diversity_data_secondary_split[[df]],ggplot2::aes(x=Window,y=gini_simpson,
                                                                                            color=Gene))+
                   ggplot2::ylim(0,1)+
                   ggplot2::theme(axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  panel.border = element_rect(color="black",fill = NA))+
                   scale_color_npg()+
                   scale_x_discrete(expand=expansion(add=50))+
                   facet_wrap(~patient_id+clinical_outcome)
           })
all_diversity_plots_secondary <-  ggarrange(plotlist = diversity_plots_secondary,nrow=5,ncol=2,common.legend = TRUE,labels = 'AUTO')

ggexport(all_diversity_plots_secondary,filename = "Output/all_diversity_plots_together_secondary.pdf",width = 8.27,height = 11.69)

#Looking at hiv and whether it is having an effect on diversity
hiv_analysis <- primary_infection_only

hiv_analysis <- hiv_analysis %>%
    mutate(`HIV status`=case_when(`HIV status`=="0" ~"negative",
                         TRUE~"positive"))
mean_of_each_patient_hiv <- hiv_analysis %>%
    group_by(patient_id,`HIV status`,clinical_outcome) %>%
    summarise(mean_by_patient=mean(gini_simpson,na.rm = TRUE))

patient_boxplot_hiv <- ggboxplot(mean_of_each_patient_hiv, y="mean_by_patient",x="HIV status",
                                       color="HIV status",outlier.shape=NA,
                                       xlab = "HIV status",ylab="Average Gini-Simpson Index",legend="none",
                                 add = "jitter",palette = c("#5C88DAFF","#CC0C00FF"))+
    stat_compare_means(label="p.format",label.x = 1.5)

mean_by_region_hiv <- hiv_analysis %>%
    group_by(Gene,`HIV status`,clinical_outcome,patient_id,) %>%
    summarise(mean_by_region=mean(gini_simpson,na.rm=TRUE))
mean_by_region_hiv$`HIV status` <-factor(mean_by_region_hiv$`HIV status`)
#Ensure genes are in correct order
mean_by_region_hiv$Gene <- factor(mean_by_region_hiv$Gene,
                                  levels = c("Core","E1","E2","P7","NS2",
                                             "NS3","NS4A","NS4B","NS5A","NS5B"))

region_plot_hiv <- ggboxplot(mean_by_region_hiv, y="mean_by_region",x="Gene",
                             color="HIV status",add.params = list(alpha = 0.3),outlier.shape = NA,
                             xlab = "Genome Regions",ylab="Average Gini-Simpson Index",legend.title="HIV status",
                             add = "jitter",palette = c("#5C88DAFF","#CC0C00FF"))+
    stat_compare_means(aes(group=`HIV status`),label="p.format",label.y = 0.7)+
    ylim(0,0.75)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
region_plot_hiv_facet <-facet(region_plot_hiv,facet.by = "Gene",scales = "free_x",nrow = 2,ncol = 5)

#Only looking at E2
e2_hiv <- mean_by_region_hiv %>%
    filter(Gene=="E2")
region_plot_hiv_e2 <- ggboxplot(e2_hiv, y="mean_by_region",x="Gene",
                             color="HIV status",add.params = list(alpha = 0.3),outlier.shape = NA,
                             xlab = "E2",ylab="Average Gini-Simpson Index",legend.title="HIV status",
                             add = "jitter",palette = c("#5C88DAFF","#CC0C00FF"))+
    stat_compare_means(aes(group=`HIV status`),label="p.format",label.y = 0.7)+
    ylim(0,0.75)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
region_plot_hiv_facet_E2 <-facet(region_plot_hiv_e2,facet.by = "clinical_outcome")


