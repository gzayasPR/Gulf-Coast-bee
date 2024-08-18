#!/bin/bash
cp ../project_env.sh .
source project_env.sh
# Now you can use $proj_dir, $my_code, $my_data, and $my_results in your script

PG_code="${my_code}/Population_Genomics/"
PG_results="${my_results}/Population_Genomics/"
PG_data="${my_data}/"
r_library="${my_softwares}/R_Packages/"

mkdir -p "${PG_code}"
mkdir -p "${PG_results}"
mkdir -p "${PG_data}"

echo "Project Directory: $proj_dir"
echo "Code Directory: $PG_code"
echo "Data Directory: $PG_data"
echo "Results Directory: $PG_results"

###Import data into PG_data###
###
env_file="${PG_code}/PG_project_env.sh"

# Write the environment variables to the file
cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="$proj_dir"
export my_code="$my_code"
export my_data="$my_data"
export my_results="$my_results"
export my_softwares="$my_softwares"
export meta_data="${meta_data}"
export r_library="${r_library}"
export ref_genome="${ref_genome}"
export PG_code="${PG_code}"
export PG_data="${PG_data}"
export PG_results="${PG_results}"

EOL

# Make the file executable
chmod +x "$env_file"

rm project_env.sh

