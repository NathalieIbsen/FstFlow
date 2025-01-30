from curses import echo
from gwf import AnonymousTarget
import os


#####################################################################
# #Templates for RUNNING Fst on all unique combination of indivuduals
# within a certain sampling setup (in workflow.py) . #
#####################################################################

def Create_all_unique(VCF: str, Temp_path: str, Populations: list, no_indv_list: list, Txt_path: str, results_path: str, REP_INFO: str, CSV_path:str):
    """Create unique sample comparisons, filter VCFs, and calculate FST."""
    # Convert no_indv_list to space-separated values for the Bash spec """ 
    no_indv_values = " ".join(map(str, no_indv_list))
    inputs = {'genome_vcf': VCF}
    outputs = {'no_indv': f"{CSV_path}complied_{Populations[0]}_{Populations[1]}_indiv_{no_indv_list[-1]}.csv"
               }
    options = {
        'cores': 2,
        'memory': '32g',
        'walltime': '24:00:00'
    }
    spec = f"""
                    # STEP 1.1 - getting variables  -----
    # Input parameters
    pop1_file="{Txt_path}pop{Populations[0]}.txt"
    pop2_file="{Txt_path}pop{Populations[1]}.txt"

    # Extract population names from file names
    pop1_id=$(basename "$pop1_file" | sed -E 's/^pop(.*)\\.txt$/\\1/')
    pop2_id=$(basename "$pop2_file" | sed -E 's/^pop(.*)\\.txt$/\\1/')

    # Read IDs from population files into arrays
    mapfile -t pop1_ids < "$pop1_file"
    mapfile -t pop2_ids < "$pop2_file"

    for no_indv in {no_indv_values}; do
        output_dir="{Temp_path}/${{pop1_id}}${{pop2_id}}/${{no_indv}}/"

        mkdir -p ${{output_dir}}  # Create Temp_path directory if it doesn't exist
        echo "Processing with no_indv = $no_indv" #process logging

        # Function to generate combinations of the size = no_indv
        generate_combinations() {{
            local size=$1
            shift
            local elements=("$@")
            local result=()
            local n=${{#elements[@]}}

            if (( size == 0 )); then
                echo ""
                return
            fi

            if (( size == 1 )); then
                for element in "${{elements[@]}}"; do
                    result+=("$element")
                done
            else
                for (( i=0; i<=n-size; i++ )); do
                    local head="${{elements[i]}}"
                    local tail_combinations
                    tail_combinations=$(generate_combinations $((size-1)) "${{elements[@]:$((i+1))}}")
                    for combination in $tail_combinations; do
                        result+=("$head,$combination")
                    done
                done
            fi

            echo "${{result[@]}}"
        }}

        ## STEP 1.2 Generate combinations for each population  -----

        pop1_combinations=($(generate_combinations $no_indv "${{pop1_ids[@]}}"))
        pop2_combinations=($(generate_combinations $no_indv "${{pop2_ids[@]}}"))

        echo "pop1_combinations: ${{pop1_combinations[@]}}"
        echo "pop2_combinations: ${{pop2_combinations[@]}}"

        # Generate output files for all combinations
        rep=1
        for comb1 in "${{pop1_combinations[@]}}"; do
            for comb2 in "${{pop2_combinations[@]}}"; do
                echo "Processing combination $comb1 and $comb2"

                # Format combinations for output
                pop1_list=$(echo "$comb1" | tr ',' '\\n')
                pop2_list=$(echo "$comb2" | tr ',' '\\n')

                # Construct output file name
                output_file="${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.txt"

                # Write combinations to file
                {{
                    echo "$pop1_list"
                    echo "$pop2_list"
                }} > "$output_file"

                ## STEP 3  -  Perform VCF filtering ----

                vcftools --vcf {VCF} --keep $output_file --out ${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}} --recode

                bcftools sort ${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.recode.vcf -Oz -o ${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.recode.vcf.gz
                bcftools index ${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.recode.vcf.gz

                vcftools --gzvcf ${{output_dir}}${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.recode.vcf.gz --max-missing 1 --out ${{output_dir}}no_miss_${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}} --recode

                #Split each pop***.txt file into their "original" division, for Fst estimation
                input_file="$output_file"
                total_lines=$(wc -l < "$input_file")
                half=$((total_lines / 2))

                if (( total_lines % 2 != 0 )); then
                    echo "ERROR: Significant Warning: The input file $input_file does not have an even number of lines. \\
                     There needs to be equal sample-size Skipping..."
                    continue
                fi

                # Construct paths for the first and bottom halves
                first_half_file="${{output_dir}}${{pop1_id}}_rep${{rep}}.txt"
                bottom_half_file="${{output_dir}}${{pop2_id}}_rep${{rep}}.txt"

                # Extract the top half of the file
                head -n "$half" "$input_file" > "$first_half_file"

                # Extract the bottom half of the file
                tail -n "$half" "$input_file" > "$bottom_half_file"

                
                # STEP 4  -  Estimate Fst ------
                vcftools --vcf ${{output_dir}}no_miss_${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.recode.vcf \\
                --weir-fst-pop "$first_half_file" --weir-fst-pop "$bottom_half_file" \\
                --out {results_path}"${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}" \\
                > ${{output_dir}}"${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.out" \\
                2> ${{output_dir}}"${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.err" \\
                
                # STEP 5 - Handle logs and results
                FST={results_path}"${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.weir.fst"

                gwf_Out=${{output_dir}}"${{pop1_id}}_${{pop2_id}}_no_indv${{no_indv}}_rep${{rep}}.err"
                # Print out two logs# 
                if [ $? -eq 0 ]; then

                    #### First log: individual FST.log#
                    
                    echo "Replicate {REP_INFO} - Replicate ${{rep}}: FST calculation completed for ${{no_indv}} individuals." > "$FST"FST.log
                    #input sample info
                    cat "$output_file" >> "$FST"FST.log

                    #Count occurrences in the 3rd column of the FST file
                    #count no of estimates by a certain quality
                    nan_count=$(awk '{{if($3 ~ /-nan/) count++}} END{{print count+0}}' "$FST")
                    gt_zero_count=$(awk '{{if($3 > 0 && $3 != "-nan") count++}} END{{print count+0}}' "$FST")
                    lt_zero_count=$(awk '{{if($3 < 0) count++}} END{{print count+0}}' "$FST")
                    #Count and log the total number of lines in the FST file
                    total_lines=$(wc -l < "$FST")

                    # Extract weighted mean from the output of vcftools (gwf_Out)
                    fst_value=$(grep "Weir and Cockerham weighted Fst estimate" "$gwf_Out" | awk -F ': ' '{{print $2}}' | awk '{{print $1}}')
            
                    # Calculate means based on the weir.fst pr. site file
                    valid_data=($(awk '$3 != "-nan" {{print $3}}' "$FST"))
                    # Mean - All negative values are set to zero
                    mean_neg_to_zero=$(echo "${{valid_data[@]}}" | tr ' ' '\\n' | awk '{{if ($1 < 0) $1=0; sum+=$1; count++}} END {{if (count > 0) print sum/count; else print "NaN"}}')
                    # Mean - All negative values > -0.10 are set to zero, any number < -0.10 are excluded from the mean
                    mean_above_threshold=$(echo "${{valid_data[@]}}" | tr ' ' '\\n' | awk '{{if ($1 >= -0.10 && $1 < 0) $1=0; if ($1 >= -0.10) {{sum+=$1; count++}}}} END {{if (count > 0) print sum/count; else print "NaN"}}')
                    
                    # Log the counts into individual FST.log (>> "$FST"FST.log) and sample info in the .log
                    echo "Total number of lines in FST file: $total_lines" >> "$FST"FST.log
                    echo "Count of -nan in 3rd column: $nan_count" >> "$FST"FST.log
                    echo "Count of values > 0 (excluding -nan) in 3rd column: $gt_zero_count" >> "$FST"FST.log
                    echo "Count of values < 0 in 3rd column: $lt_zero_count" >> "$FST"FST.log
                    echo "______Mean estimates Replicate replicate ${{rep}} {Populations[0]}_{Populations[1]}__________" >> "$FST"FST.log
                    echo "Weir and Cockerham weighted Fst estimate: $fst_value" >> "$FST"FST.log
                    echo "Mean - All negative valoues are set to zero: $mean_neg_to_zero" >> "$FST"FST.log
                    echo "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean: $mean_above_threshold" >> "$FST"FST.log
                else
                    echo "ERROR: FST calculation failed for ${{no_indv}} individuals." >> "$FST"FST.log
                    cat "$output_file" >> "$FST"FST.log
                fi 
                
                # Increment replicate count
                rep=$((rep + 1))
            done
        done
    done
                #STEP 2 joining the replicates 

    find {results_path} -type f -exec chmod 777 {{}} +

    echo "permission to {results_path}"

    for no_indv in {no_indv_values}; do

        # Define output file for this specific number of individuals
        OUTPUT_FILE="{CSV_path}complied_{Populations[0]}_{Populations[1]}_indiv_${{no_indv}}.csv"
            # Print the header for the CSV file
            echo "Replicate_no,Weir_and_Cockerham_Fst,Mean_Negatives_Set_to_Zero,Mean_Threshold_Exclusion,no_sites,no_monomorph" > "$OUTPUT_FILE"
            # Loop over all .fst log files matching the pattern
            for FILE in {results_path}{Populations[0]}_{Populations[1]}_no_indv${{no_indv}}_rep*.weir.fstFST.log; do
                #debug:
                echo "Processing: $FILE" 
                    # Check if the file exists
                if [ ! -f "$FILE" ]; then
                    echo "File does not exist or is inaccessible: $FILE"
                    continue
                fi

                # Extract replicate number from filename (assuming pattern 'rep_int')
                rep=$(basename "$FILE" | grep -oP 'rep\\K\\d+')
                echo "Extracted replicate number: $rep"

                # Extract FST statistics from the file
                Weir_Fst=$(grep "Weir and Cockerham weighted Fst estimate" "$FILE" | awk -F': ' '{{print $2}}')
                Mean_Set_to_Zero=$(grep "Mean - All negative valoues are set to zero" "$FILE" | awk -F': ' '{{print $2}}')
                Mean_Threshold_Exclusion=$(grep "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean" "$FILE" | awk -F': ' '{{print $2}}')
                no_sites=$(grep "Total number of lines in FST file" "$FILE" | awk -F': ' '{{print $2}}')
                no_monomorphic=$(grep "Count of -nan in 3rd column" "$FILE" | awk -F': ' '{{print $2}}')
                # Append the extracted data to the output CSV
                #Debug
                #echo "$rep,$Weir_Fst,$Mean_Set_to_Zero,$Mean_Threshold_Exclusion,$no_sites,$no_monomorphic"
                echo "$rep,$Weir_Fst,$Mean_Set_to_Zero,$Mean_Threshold_Exclusion,$no_sites,$no_monomorphic" >> "$OUTPUT_FILE"
            done

            # Sort the CSV file by Replicate_no (first column), excluding the header
            awk 'BEGIN{{FS = OFS = ","}} {{if (NR==1) {{print $0; next}} print $0 | "sort -t, -k1,1n "}}' "$OUTPUT_FILE" > "$OUTPUT_FILE.sorted" && mv "$OUTPUT_FILE.sorted" "$OUTPUT_FILE"
    done  
    echo "Data extraction completed. Check the output in the output files within {CSV_path}."  
    
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def Get_FST_compiled_CSV(CSVs: list, CSV_path: str):
    """Extract vcftools weighted mean, max and min of each replication run and calculate the mean we are interested in."""

    inputs = {
        'FST_CSV': CSVs
    }
    outputs = {
        'FST_csv': f"{CSV_path}/concatenated_runs_all_populations_all_replications.csv"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    # Define the shell script to be executed
    spec = """
    #############CONCAT all csv############
    # Define the directory where the CSV files are located 
    INPUT_DIR={CSV_path}
    OUTPUT_FILE="{CSV_path}/concatenated_runs_all_populations_all_replications.csv"

    # Initialize a flag to check if it's the first file
    first_file=true
    # Loop through each CSV file in the directory
    for file in "$INPUT_DIR"/complied_*.csv; do
        # Extract the population identifier from the filename
        # Example: complied_SPA_HOR_indiv_2.csv -> SPA_HOR
        population=$(basename "$file" | awk -F'_' '{{print $2 "_" $3}}' | sed 's/.csv//')

        # Extract the no_of_indv from the filename (the last number before .csv)
        no_of_indv=$(basename "$file" | grep -o '[0-9]\+' | tail -1)

        # Print header for the output file using the header from the first file
        if [ "$first_file" = true ]; then
            header=$(head -n 1 "$file")  # Capture the full header directly
            echo "population,no_of_indv,$header" > "$OUTPUT_FILE"
            first_file=false  # Set the flag to false after the first file is processed
        fi
        # Read the CSV file and append data to the output file
        # Skip the header (first line) of each file, and add the population identifier and no_of_indv as the first two columns
        awk -v pop="$population" -v indiv="$no_of_indv" 'NR > 1 {{print pop "," indiv "," $0}}' "$file" >> "$OUTPUT_FILE"
    done
    echo "Concatenation completed. Output saved in $OUTPUT_FILE."

    """.format(
        CSV_path=CSV_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def Clean(CSV: str, Temp_path: str, Out_path: str, CSV_path: str):
    """Cleaning up post run, remove all temperoary files"""
    inputs = {
        'CSV': CSV
    }
    outputs = { 'FST_csv': f"{CSV_path}/concatenated_runs_all_populations_all_replications.csv"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    # Define the shell script to be executed
    spec = """
    echo "Removing {Temps_path}..."
    rm -rf {Temps_path}
    echo "Removing {Out_path}..."
    rm -rf {Out_path}
    """.format(
        Temps_path=Temp_path,
        Out_path=Out_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#find {Out_path} -type f -name "*.fst" -exec rm -f {} \;