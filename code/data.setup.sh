#!/bin/bash

# Source the environment variables
source project_env.sh

# Define the paths to the reference genome and metadata
reference_genome=$1
metadata=$2

# Create the necessary directories in the data folder
mkdir -p $my_data/reference_genome/

# Copy the reference genome to the data folder
cp ${reference_genome} $my_data/reference_genome/.

# Get the basename of the copied reference genome
ref_genome="$my_data/reference_genome/$(basename $reference_genome)"

# Copy the metadata file to the data folder
cp ${metadata} $my_data/

# Define the full path of the copied metadata file
meta_data="$my_data/$(basename ${metadata})"

# Print the full paths of the copied files
echo "$ref_genome"
echo "$meta_data"

# Append the ref_genome and meta_data variables to project_env.sh
echo "export ref_genome=\"$ref_genome\"" >> "$my_code/project_env.sh"
echo "export meta_data=\"$meta_data\"" >> "$my_code/project_env.sh"


