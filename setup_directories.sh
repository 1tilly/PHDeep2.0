#!/bin/bash

# Base directories
BASE_DIRS=("data/bed_files" "data/vcf_files" "data/reference_genome" "src/data_loading" "src/data_processing" "src/models/base_model" "src/models/model2" "src/models/bert_models/bert_model1" "src/prediction" "src/post_prediction" "src/statistical_testing" "src/utils" "src/experiments" "config" "results")

# Files to be created within specific directories
declare -A FILES
FILES["src/data_loading"]="bed_loader.py vcf_loader.py genome_loader.py gff_loader.py"
FILES["src/data_processing"]="bed_to_training.py vcf_processing.py"
FILES["src/models/base_model"]="architecture.py train.py"
FILES["src/models/cnn_model1"]="architecture.py train.py"
FILES["src/models/bert_models/bert_model1"]="architecture.py train.py tokeniser.py vocabulary.txt"
FILES["src/prediction"]="predict.py"
FILES["src/post_prediction"]="aggregation.py"
FILES["src/statistical_testing"]="skat_o_test.py"
FILES["src/utils"]="plotting.py"
FILES["src/experiments"]="base_model_config.yaml cnn_model1_config.yaml bert_model1_config.yaml"
FILES["config"]="paths.py"
FILES["."]="README.md requirements.txt main.py"

# Create base directories
for dir in "${BASE_DIRS[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
  fi
done

# Create specific files within directories
for dir in "${!FILES[@]}"; do
  IFS=' ' read -ra ADDR <<< "${FILES[$dir]}"
  for file in "${ADDR[@]}"; do
    if [ ! -f "$dir/$file" ]; then
      touch "$dir/$file"
    fi
  done
done

echo "Directory and file structure set up completed."

