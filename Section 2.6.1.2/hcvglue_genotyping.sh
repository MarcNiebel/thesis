#!/bin/bash

#Created by Marc Niebel Feb 2021
#This is using the HCV-GLUE framework to genotype based on a file containing consensus sequences
#Output is the genotype and subtype

gluetools.sh -p cmd-result-format:tab -EC --inline-cmd project hcv module maxLikelihoodGenotyper genotype file -f $1 > all_genotyping-$$
