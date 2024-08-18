#!/bin/bash


reference_genome=/project/90daydata/beenome100/hesperapis_oraria_genomics/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
metadata=/project/90daydata/beenome100/hesperapis_oraria_genomics/Samples.Metadata.csv

proj_dir="$(pwd)"
my_code="$(pwd)/code"
my_data="$(pwd)/data"
my_results="$(pwd)/results"
my_softwares="$(pwd)/softwares"

mkdir -p $my_code
mkdir -p $my_data
mkdir -p $my_results
mkdir -p $my_softwares

env_file="$my_code/project_env.sh"

cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="$(pwd)"
export my_code="$(pwd)/code"
export my_data="$(pwd)/data"
export my_results="$(pwd)/results"
export my_softwares="$(pwd)/softwares"
EOL

# Make the file executable
chmod +x "$env_file"


cd $my_code

echo "Starting refgenome setup"
bash data.setup.sh ${reference_genome} ${metadata}
echo "Fineshed refgenome setup"

echo "Starting software Setup"
bash Software.setup.sh


bash code.dir.setup.sh