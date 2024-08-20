#!/bin/bash
cp ../project_env.sh .
source project_env.sh
# Now you can use $proj_dir, $my_code, $my_data, and $my_results in your script

ANGSD_code="${my_code}/ANGSD_code/"
ANGSD_results="${my_results}/ANGSD/"
ANGSD_data="${my_data}/"
r_library="${my_softwares}/R_Packages/"

mkdir -p "${ANGSD_code}"
mkdir -p "${ANGSD_results}"
mkdir -p "${ANGSD_data}"

echo "Project Directory: $proj_dir"
echo "Code Directory: $ANGSD_code"
echo "Data Directory: $ANGSD_data"
echo "Results Directory: $ANGSD_results"

###Import data into ANGSD_data###
###
env_file="${ANGSD_code}/ANGSD_project_env.sh"

# Write the environment variables to the file
cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="$proj_dir"
export my_code="$my_code"
export my_data="$my_data"
export my_results="$my_results"
export my_softwares="$my_softwares"
export my_docs="$my_docs"
export meta_data="${meta_data}"
export r_library="${r_library}"
export ref_genome="${ref_genome}"
export ANGSD_code="${ANGSD_code}"
export ANGSD_data="${ANGSD_data}"
export ANGSD_results="${ANGSD_results}"

EOL

# Make the file executable
chmod +x "$env_file"

rm project_env.sh

