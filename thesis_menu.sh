#!/bin/bash

COLORBLUE=$'\e[1;44m'
COLORYELLOW=$'\e[1;43m'
COLORRED=$'\e[1;41m'
NCOLOR='\033[0m'

## SCRIPTS

thesis_rscript="thesis_analysis.R"
thesis_pyscript="thesis_analysis.py"

for study in configs/*; do
  for conf in $study/*; do

    # GET PARAMETERS FOR CONFIG
    . $conf

    if [[ ${condition} == 1 ]]; then

      eval "$(conda shell.bash hook)"

      printf "${COLORBLUE}ACTIVATING CONDA ENV${NCOLOR}\n"

      conda activate thesis_r

      printf "${COLORBLUE}CREATING STUDY FOLDER IN results/${NCOLOR}\n"

      mkdir -p results/${study_name}

      printf "${COLORBLUE}STARTING RSCRIPT${NCOLOR}\n"
      Rscript $thesis_rscript --study_name "${study_name}" \
                              --biom_file "${study_biom}" \
                              --mapping_file "${mapping_file}" \
                              --taxonomy_file "${taxonomy_file}" \
                              --tree_file "${tree_file}" \
                              --class_variable "${class_variable}" \
                              --classes "${classes}" \
                              --taxa_rank "${taxa_rank}" \
                              --transformation "${transformation}" \
                              --distance_measure "${distance_measure}" \
                              --phase1 TRUE \
                              --phase2 FALSE

      printf "${COLORBLUE}PHASE 1 RESULTS DONE${NCOLOR}\n"

      printf "${COLORYELLOW}ACTIVATING CONDA ENV${NCOLOR}\n"
      conda activate thesis_ml

      printf "${COLORYELLOW}STARTING PYSCRIPT${NCOLOR}\n"
      python3 $thesis_pyscript --study_name "${study_name}" \
                               --taxa_rank "${taxa_rank}" \
                               --transformation "${transformation}" \
                               --distance_measure "${distance_measure}" \

      conda activate thesis_r

      printf "${COLORYELLOW}STARTING RSCRIPT${NCOLOR}\n"
      Rscript $thesis_rscript --study_name "${study_name}" \
                              --biom_file "${study_biom}" \
                              --mapping_file "${mapping_file}" \
                              --taxonomy_file "${taxonomy_file}" \
                              --tree_file "${tree_file}" \
                              --class_variable "${class_variable}" \
                              --classes "${classes}" \
                              --taxa_rank "${taxa_rank}" \
                              --transformation "${transformation}" \
                              --distance_measure "${distance_measure}" \
                              --phase1 FALSE \
                              --phase2 TRUE

    elif [[ ${condition} == 2 ]]; then

      eval "$(conda shell.bash hook)"
      printf "${COLORYELLOW}ACTIVATING CONDA ENV${NCOLOR}\n"
      conda activate thesis_ml

      printf "${COLORYELLOW}STARTING PYSCRIPT${NCOLOR}\n"
      python3 $thesis_pyscript --study_name "${study_name}" \
                               --taxa_rank "${taxa_rank}" \
                               --transformation "${transformation}" \
                               --distance_measure "${distance_measure}" \

      conda activate thesis_r

      printf "${COLORYELLOW}STARTING RSCRIPT${NCOLOR}\n"
      Rscript $thesis_rscript --study_name "${study_name}" \
                              --biom_file "${study_biom}" \
                              --mapping_file "${mapping_file}" \
                              --taxonomy_file "${taxonomy_file}" \
                              --tree_file "${tree_file}" \
                              --class_variable "${class_variable}" \
                              --classes "${classes}" \
                              --taxa_rank "${taxa_rank}" \
                              --transformation "${transformation}" \
                              --distance_measure "${distance_measure}" \
                              --phase1 FALSE \
                              --phase2 TRUE

    else

      printf "${COLORRED}PLEASE CHECK 'thesis_parameters.config' CONFIG FILE FOR CONDITION SELECTED${NCOLOR}\n"

    fi

 done
done
