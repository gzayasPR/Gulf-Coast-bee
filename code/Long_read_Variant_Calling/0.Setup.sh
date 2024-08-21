#!/bin/bash
cp ../project_env.sh .
source project_env.sh
# Now you can use $proj_dir, $my_code, $my_data, and $my_results in your script

LR_VC_code="${my_code}/Long_read_Variant_Calling/"
LR_VC_results="${my_results}/Long_read_Variant_Calling/"
LR_VC_data="${my_data}/"
r_library="${my_softwares}/R_Packages/"

mkdir -p "${LR_VC_code}"
mkdir -p "${LR_VC_results}"
mkdir -p "${LR_VC_data}"

echo "Project Directory: $proj_dir"
echo "Code Directory: $LR_VC_code"
echo "Data Directory: $LR_VC_data"
echo "Results Directory: $LR_VC_results"

###Import data into LR_VC_data###
###
env_file="${LR_VC_code}/LR_VC_project_env.sh"

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
export LR_VC_code="${LR_VC_code}"
export LR_VC_data="${LR_VC_data}"
export LR_VC_results="${LR_VC_results}"
EOL

# Make the file executable
chmod +x "$env_file"

rm project_env.sh
