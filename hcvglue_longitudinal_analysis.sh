#!/bin/bash

#Created by Marc Niebel November 2020
#This script uses the HCV-GLUE framework to investigate the amino acid dynamics for longitudinally sampled patients

#Supply a list of folder ids where bams are stored

#For loop which reads in each patient id indidvidually

cat $1 | while read patient_id || [[ -n $patient_id ]];
do
        gluetools.sh --inline-cmd project hcv module thacBamUtilities invoke-function aminoAcidDynamicsToFile $patient_id
done
