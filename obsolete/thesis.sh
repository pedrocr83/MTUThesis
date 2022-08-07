#!/bin/bash


kraken2="bin/kraken2/./kraken2" # kraken2 software
bracken="bin/Bracken-master/./bracken" # bracken software
guthealth="bin/kraken2/guthealth/" # kraken2 database to be used
thesis_phase1="thesis_phase1.R" # R script to generate metrics


################################
#      ACTIVATE CONDA ENV      #
################################

eval "$(conda shell.bash hook)"
echo "ACTIVATING CONDA ENV"
conda activate thesis

printf "KRAKEN CONFIDENCE THRESHOLD? (0-1)\n"
read threshold_num

###############################
#           PHASE 1           #
###############################

printf "\nCREATING PHASE 1 METRICS\n"

cd src

Rscript $thesis_phase1 -s "${sample_id}" -b "../${i%%.*}.biom"
