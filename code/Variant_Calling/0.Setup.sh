#!/bin/bash
cp ../project_env.sh .
source project_env.sh
# Now you can use $proj_dir, $my_code, $my_data, and $my_results in your script

VC_code="${my_code}/Variant_Calling/"
VC_results="${my_results}/Variant_Calling/"
VC_data="${my_data}/"
r_library="${my_softwares}/R_Packages/"

mkdir -p "${VC_code}"
mkdir -p "${VC_results}"
mkdir -p "${VC_data}"

echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

###Import data into VC_data###
###
env_file="${VC_code}/VC_project_env.sh"

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
export VC_code="${VC_code}"
export VC_data="${VC_data}"
export VC_results="${VC_results}"
EOL

# Make the file executable
chmod +x "$env_file"

rm project_env.sh
