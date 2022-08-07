#!/bin/bash

COLORBLUE=$'\e[1;44m'
COLORYELLOW=$'\e[1;43m'
COLORRED=$'\e[1;41m'
NCOLOR='\033[0m'

## SCRIPTS

thesis_rscript="thesis_analysis.R"
thesis_pyscript="thesis_analysis.py"

# GET PARAMETERS FOR CONFIG
. thesis_parameters.config

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
                           --knn_nneighbours "${knn_nneighbours}" \
                           --outlier_value "${outlier_value}"

  conda activate thesis_r

  printf "${COLORYELLOW}STARTING RSCRIPT${NCOLOR}\n"
  Rscript $thesis_rscript --study_name "${study_name}" \
                          --biom_file "${study_biom}" \
                          --mapping_file "${mapping_file}" \
                          --taxonomy_file "${gg_taxonomy}" \
                          --tree_file "${gg_tree}" \
                          --class_variable "${class_variable}" \
                          --classes "${classes}" \
                          --taxa_rank "${taxa_rank}" \
                          --transformation "${transformation}" \
                          --distance_measure "${distance_measure}" \
                          --distance_measure_results "${distance_measure_results}" \
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
                           --knn_nneighbours "${knn_nneighbours}" \
                           --outlier_value "${outlier_value}"

  conda activate thesis_r

  printf "${COLORYELLOW}STARTING RSCRIPT${NCOLOR}\n"
  Rscript $thesis_rscript --study_name "${study_name}" \
                          --biom_file "${study_biom}" \
                          --mapping_file "${mapping_file}" \
                          --taxonomy_file "${gg_taxonomy}" \
                          --tree_file "${gg_tree}" \
                          --class_variable "${class_variable}" \
                          --classes "${classes}" \
                          --taxa_rank "${taxa_rank}" \
                          --transformation "${transformation}" \
                          --distance_measure "${distance_measure}" \
                          --distance_measure_results "${distance_measure_results}" \
                          --phase1 FALSE \
                          --phase2 TRUE

else

  printf "${COLORRED}PLEASE CHECK 'thesis_parameters.config' CONFIG FILE FOR CONDITION SELECTED${NCOLOR}\n"

fi



