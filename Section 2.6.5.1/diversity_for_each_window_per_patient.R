#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

#Created in November 2021
#Purpose is  to calculate the diversity measurement(Gini-Simpson Index)
#for all windows across the genome for an individual patient which has been written
#to be used on a folder of multiple sequence alignments

#libraries needed
suppressPackageStartupMessages({
    library(stringi)
    library(QSutils)
    library(gtools)
    library(dplyr)
})
alignment_files <- list.files(args[1],pattern = "*.goodfna",
                              full.names = TRUE)
#Patient id extracted
patient_id <- stri_extract(alignment_files[1], regex="\\w\\d{1,}")

#Making a list of sequence alignment files
sequence_files <- lapply(alignment_files,readDNAStringSet)    

#Gini-Simpson Index diversity measurements####
#Generating the haplotypes for each alignment
reads_collapsed <- lapply(sequence_files,Collapse)
#Pulling out the frequencies
haplotype_freq <- sapply(reads_collapsed,"[",1)
gini_simpson <-sapply(haplotype_freq,GiniSimpson)
gini_simpson_df <- as.data.frame(gini_simpson)

#Adding on the amino acid coordinates to the diversity dataframe####
#Extracting the amino acid coordinates from each file and putting them in a list
amino_acid_cooridnates <- lapply(alignment_files, function(name) 
    as.character(stri_match_last(name, regex="\\d{1,}_\\d{1,}")))
#Needed as a character vector
amino_acid_cooridnates <- as.character(amino_acid_cooridnates)
#Added on the diversity measurement dataframe
gini_simpson_df$amino_acid_windows <- amino_acid_cooridnates
#Ordered them based on window coordinates
gini_simpson_df_ordered <- gini_simpson_df[mixedorder(gini_simpson_df$amino_acid_windows),]

#Adding on which gene window originates(H77)from and also have missing data available####
#Reading in relevant file
window_coordinates_associated_genes <- read.csv("Data/Window_coordinates_nine_amino_acid_H77.csv")
#Adding on diversity measurements ensuring that missing data is maintained
adding_genes_diversity <- left_join(window_coordinates_associated_genes,
                                    gini_simpson_df_ordered,by=c("Window" ="amino_acid_windows"))
#Adding patient id
adding_genes_diversity$patient_id <- patient_id

write.csv(adding_genes_diversity,paste0("Output/Baseline diversity/",patient_id,"_all_window_diversity.csv"),
          row.names = FALSE)
