#!/bin/bash

## FILES

dataset_url="http://metagenome.cs.umn.edu/public/MLRepo/datasets.tar.gz" # MLRepo Datasets found on https://knights-lab.github.io/MLRepo/
thesis_env_yml="thesis.yml"
taxonomy_url="https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_taxonomy.txt.gz"
tree_url="https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_otus_99_annotated.tree.gz"

dataset_tar="datasets.tar.gz"
taxonomy_gz="gg_13_5_taxonomy.txt.gz"
tree_gz="gg_13_5_otus_99_annotated.tree.gz"

## SCRIPTS

preparation_R="preparation.R"

## VARIABLES

thesis_env="thesis2"
taxonomy_file="gg_13_5_taxonomy.txt"
tree_file="gg_13_5_otus_99_annotated.tree"



download_and_install() {
#  echo "DOWNLOADING DATASETS AND RELEVANT FILES"
#  echo ""
#  wget ${dataset_url}
#  wget ${taxonomy_url}
#  wget ${tree_url}
#  tar xvzf ${dataset_tar}
#  gunzip ${taxonomy_gz}
#  gunzip ${tree_gz}
  
#  echo ""
#  echo "INSTALLING MAMBA PACKAGE (FOR CONDA SPEED)"
#  echo ""
#  conda install -y mamba -n base -c conda-forge
#  echo "CREATING THESIS CONDA ENVIRONMENT"

  eval "$(conda shell.bash hook)"
  
  #conda env create --quiet -f ${thesis_env_yml}
  echo "ACTIVATING THESIS CONDA ENVIRONMENT"
  echo ""
  conda activate ${thesis_env}
  echo "INSTALLING R ADDITIONAL PACKAGES"
  echo ""
  Rscript $preparation_R
  
}

run_thesis_pipeline() {
  
  printf "CHOOSE DATASET?\n"
  read dataset_folder
  
  printf "CHOOSE DATASET?\n"
  read dataset_folder

  printf "CHOOSE DATASET?\n"
  read dataset_folder
  
  printf "CHOOSE DATASET?\n"
  read dataset_folder
  
}

press_enter() {
  echo ""
  echo -n "	Press Enter to continue "
  read
  clear
}

incorrect_selection() {
  echo "Incorrect selection! Try again."
}

until [ "$selection" = "0" ]; do
  clear
  echo ""
  echo "    	1  -  Download Data, create Conda Environment and Install Packages"
  echo "    	2  -  Run Thesis Experiment Pipeline"
  echo "    	0  -  Exit"
  echo ""
  echo -n "  Enter selection: "
  read selection
  echo ""
  case $selection in
    1 ) clear ; download_and_install ; press_enter ;;
    2 ) clear ; run_thesis_pipeline ; press_enter ;;
    0 ) clear ; exit ;;
    * ) clear ; incorrect_selection ; press_enter ;;
  esac
done
