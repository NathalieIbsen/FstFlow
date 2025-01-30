from gwf import AnonymousTarget
import os

#####################################################################
# #Templates for RUNNING Fst on a random sampled no of indivuduals. #
#####################################################################


def Random_draw_indv(VCF: str, Temp_path: str, Populations: list, no_indv: int, Txt_path: str):
    """Draw a defined number of individuals from VCF from two populations.
    2)Sort, index and write a popualtion file, 
    3)Filter out missing positions,
    4) Merge population -> vcf, & write .txt"""
    inputs = {'genome_vcf': VCF}
    outputs = {
        'Temp_vcf_1': f"{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz",
        'Temp_vcf_0': f"{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz",
        'Indvi_txt_1': f"{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt",
        'Indvi_txt_0': f"{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt",
        'Merged_vcf': f"{Temp_path}/no_miss_{Populations[0]}_{Populations[1]}_Merged_{no_indv}.recode.vcf"
    }
    options = {
        'cores': 6,
        'memory': '16g',
        'walltime': '04:00:00'
    }

    spec = f"""
    mkdir -p {Temp_path}  # Create directory, if it exists - skip

    for pop in {' '.join(Populations)}; do
        vcftools --vcf {VCF} --keep {Txt_path}/pop${{pop}}.txt \\
        --max-indv {no_indv} \\
        --out {Temp_path}/${{pop}}_Tempory_{no_indv} --recode

        bcftools sort {Temp_path}/${{pop}}_Tempory_{no_indv}.recode.vcf -Ob \\
        -o {Temp_path}/${{pop}}_Tempory_{no_indv}.recode.bcf.gz

        bcftools index {Temp_path}/${{pop}}_Tempory_{no_indv}.recode.bcf.gz

        bcftools query -l {Temp_path}/${{pop}}_Tempory_{no_indv}.recode.bcf.gz \\
        > {Temp_path}/${{pop}}_Tempory_{no_indv}.txt
    done

    # Rename files to include both population names
    mv "{Temp_path}/{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz" \\
       "{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz"
    mv "{Temp_path}/{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz" \\
       "{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz"
    mv "{Temp_path}/{Populations[1]}_Tempory_{no_indv}.txt" \\
       "{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt"
    mv "{Temp_path}/{Populations[0]}_Tempory_{no_indv}.txt" \\
       "{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt"
    mv "{Temp_path}/{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz.csi" \\
       "{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz.csi"
    mv "{Temp_path}/{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz.csi" \\
       "{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz.csi"

    bcftools merge \\
        {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz \\
        {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz \\
        --threads 6 \\
        --output-type v \\
        --output {Temp_path}/{Populations[0]}_{Populations[1]}_Merged_{no_indv}.vcf

    vcftools --vcf {Temp_path}/{Populations[0]}_{Populations[1]}_Merged_{no_indv}.vcf \\
        --max-missing 1 \\
        --out {Temp_path}/no_miss_{Populations[0]}_{Populations[1]}_Merged_{no_indv} --recode
    
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


####NOT needed anymore as the above job now handles the merging. 
# 2nd Job: Merge subsampled individuals from both populations
def Merge_indv_files(Temp_vcf_1: str, Temp_vcf_0: str, Temp_path: str, Populations: list, no_indv: int):
    """Merge VCF files of subsampled individuals from two populations. + Filter missing sites"""
    inputs = {
        'Indvi_txt_1': f"{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt",
        'Indvi_txt_0': f"{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt"
    }
    outputs = {
        'Merged_bcf': f"{Temp_path}/old_no_miss_{Populations[0]}_{Populations[1]}_Merged_{no_indv}.recode.vcf"
    }
    options = {
        'cores': 3,
        'memory': '16g',
        'walltime': '04:00:00'
    }
    spec = """
    bcftools merge \\
        {Temp_vcf_1} \\
        {Temp_vcf_0} \\
        --threads 3 \\
        --output-type v \\
        --output {Temp_path}/{Populations[0]}_{Populations[1]}_Merged_{no_indv}.vcf

    vcftools --vcf {Temp_path}/{Populations[0]}_{Populations[1]}_Merged_{no_indv}.vcf \\
        --max-missing 1 \\
        --out {Temp_path}/no_miss_{Populations[0]}_{Populations[1]}_Merged_{no_indv} --recode
    
    """.format(Temp_path=Temp_path, Populations=Populations, no_indv=no_indv, Temp_vcf_0=Temp_vcf_0, Temp_vcf_1=Temp_vcf_1)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



# 3rd Job: Compute FST for the merged files
def Calc_FST_files(Merged_vcf: str, Temp_path: str, Output_path: str, scripts_path:str, Populations: list, no_indv: int, Replicate_no: int, REP_INFO: str):
    """Compute FST for subsampled individuals from two populations + calculate means + Print logs"""
    inputs = {
        'Merged_vcf': Merged_vcf
    }
    outputs = {
        'FST_file': f"{Output_path}/{REP_INFO}_ReP_{Replicate_no}_{Populations[0]}_{Populations[1]}_FST_{no_indv}.weir.fst", 
        'log': f"{Output_path}/{REP_INFO}_ReP_{Replicate_no}_{Populations[0]}_{Populations[1]}_FST_{no_indv}.weir.fstFST.log"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    spec = """
    #Run Fst software 
    vcftools --gzvcf {Merged_vcf} \\
        --weir-fst-pop {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt \\
        --weir-fst-pop {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt \\
        --out {Output_path}/{REP_INFO}_ReP_{Replicate_no}_{Populations[0]}_{Populations[1]}_FST_{no_indv}

    FST={Output_path}/{REP_INFO}_ReP_{Replicate_no}_{Populations[0]}_{Populations[1]}_FST_{no_indv}.weir.fst

    gwf_Out="{scripts_path}/.gwf/logs/FST_{Replicate_no}{Populations[0]}{Populations[1]}{no_indv}{REP_INFO}.stderr"
    # Print out two logs# 
    if [ $? -eq 0 ]; then

        #### First log: individual FST.log#
        
        echo "Replicate {REP_INFO} - Replicate {Replicate_no}: FST calculation completed for {no_indv} individuals." >> "$FST"FST.log
        #input sample info
        cat {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt >> "$FST"FST.log
        cat {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt >> "$FST"FST.log

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
        echo "______Mean estimates Replicate {Replicate_no} {Populations[0]}_{Populations[1]}__________" >> "$FST"FST.log
        echo "Weir and Cockerham weighted Fst estimate: $fst_value" >> "$FST"FST.log
        echo "Mean - All negative valoues are set to zero: $mean_neg_to_zero" >> "$FST"FST.log
        echo "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean: $mean_above_threshold" >> "$FST"FST.log

        ##Secondary run.log.. Contains info for all runes in the replication run 
        # Also log a summary in the main FST_RUN.log
        echo "Replicate {REP_INFO} - Replicate {Replicate_no}: FST calculation completed for {no_indv} individuals." >> {Output_path}{REP_INFO}FST_RUN.log
        # Second log: Append content of temporary files to FST_RUN.log
        cat {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt >> {Output_path}{REP_INFO}FST_RUN.log
        cat {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt >> {Output_path}{REP_INFO}FST_RUN.log
    else
        echo "ERROR: FST calculation failed for {no_indv} individuals." >> "$FST"FST.log
        echo "{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt" >> "$FST"FST.log
        echo "{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt" >> "$FST"FST.log

        # Error handling: Log errors and file paths if FST calculation fails
        echo "ERROR: FST calculation failed for {no_indv} individuals." >> {Output_path}FST_RUN.log
        echo "{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt" >> {Output_path}FST_RUN.log
        echo "{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt" >> {Output_path}FST_RUN.log
    fi 
    """.format(
        Merged_vcf=Merged_vcf,
        Temp_path=Temp_path, 
        Output_path=Output_path, 
        Populations=Populations, 
        no_indv=no_indv, 
        REP_INFO=REP_INFO, 
        Replicate_no=Replicate_no,
        scripts_path=scripts_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def Get_FST_means(Output_path: str, FSTs: list, no_indv_list: list, Populations: list, Input_path: str, means_path: str):
    """Extract vcftools weighted mean, max and min of each replication run and calculate the mean we are interested in."""

    inputs = {
        'FST_logs': FSTs
    }
    outputs = {
        'Fst_csv': f"{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_{no_indv_list[0]}.csv",
        'Means_file':f"{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '00:10:00'
    }
    
    # Define the shell script to be executed
    spec = """
        # Set directory containing FST files
        FST_DIR="{Input_path}"
        
        # Loop over each number in the no_indv_list
        for no_indv in {no_indv_list}; do
            # Define output file for this specific number of individuals
            OUTPUT_FILE="{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_$no_indv.csv"

            # Print the header for the CSV file
            echo "Replicate_no,Weir_and_Cockerham_Fst,Mean_Negatives_Set_to_Zero,Mean_Threshold_Exclusion,no_sites,no_monomorphic" > "$OUTPUT_FILE"
            
            # Loop over all .fst log files matching the pattern
            for FILE in "$FST_DIR"/*_{Populations[0]}_{Populations[1]}_FST_"$no_indv".weir.fstFST.log; do
                # Extract replicate number from filename (assuming pattern 'ReP_<Replicate_no>')
                Replicate_no=$(basename "$FILE" | grep -oP 'ReP_\\K\\d+')
                
                # Extract FST statistics from the file
                Weir_Fst=$(grep "Weir and Cockerham weighted Fst estimate" "$FILE" | awk -F': ' '{{print $2}}')
                Mean_Set_to_Zero=$(grep "Mean - All negative valoues are set to zero" "$FILE" | awk -F': ' '{{print $2}}')
                Mean_Threshold_Exclusion=$(grep "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean" "$FILE" | awk -F': ' '{{print $2}}')
                no_sites=$(grep "Total number of lines in FST file" "$FILE" | awk -F': ' '{{print $2}}')
                no_monomorphic=$(grep "Count of -nan in 3rd column" "$FILE" | awk -F': ' '{{print $2}}')
                # Append the extracted data to the output CSV
                echo "$Replicate_no,$Weir_Fst,$Mean_Set_to_Zero,$Mean_Threshold_Exclusion,$no_sites,$no_monomorphic" >> "$OUTPUT_FILE"
            done

            # Sort the CSV file by Replicate_no (first column), excluding the header
            awk 'BEGIN{{FS = OFS = ","}} {{if (NR==1) {{print $0; next}} print $0 | "sort -t, -k1,1n "}}' "$OUTPUT_FILE" > "$OUTPUT_FILE.sorted" && mv "$OUTPUT_FILE.sorted" "$OUTPUT_FILE"
        done  
        echo "Data extraction completed. Check the output in the output files within {Output_path}."

        # Print the data in the CSV file
        echo "no_of_indv,Mean_Weir_and_Cockerham_Fst,Min_Weir_and_Cockerham_Fst,Max_Weir_and_Cockerham_Fst,Mean_Mean_Negatives_Set_to_Zero,Min_Mean_Negatives_Set_to_Zero,Max_Mean_Negatives_Set_to_Zero,Mean_Mean_Threshold_Exclusion,Min_Mean_Threshold_Exclusion,Max_Mean_Threshold_Exclusion" > "{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
        Means_FILE="{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
        # Loop over each number in the no_indv_list
        for no_indv in {no_indv_list}; do
            CSV_FILE="{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_$no_indv.csv"
            #Make sure premissions are given
            chmod 777 "$CSV_FILE"
            # Calculate the mean, min, and max of the replicates - Store in variables 
            read Mean_Weir_and_Cockerham_Fst Min_Weir_and_Cockerham_Fst Max_Weir_and_Cockerham_Fst Mean_Mean_Negatives_Set_to_Zero Min_Mean_Negatives_Set_to_Zero Max_Mean_Negatives_Set_to_Zero Mean_Mean_Threshold_Exclusion Min_Mean_Threshold_Exclusion Max_Mean_Threshold_Exclusion < <(
                awk -F',' 'NR > 1 {{ # skip header 
                    sum2 += $2; sum3 += $3; sum4 += $4;
                    if ($2 < min2 || NR == 2) min2 = $2;
                    if ($2 > max2 || NR == 2) max2 = $2;
                    if ($3 < min3 || NR == 2) min3 = $3;
                    if ($3 > max3 || NR == 2) max3 = $3;
                    if ($4 < min4 || NR == 2) min4 = $4;
                    if ($4 > max4 || NR == 2) max4 = $4;
                    count += 1
                }} END {{
                    if (count > 0) {{
                        print (sum2 / count), min2, max2, (sum3 / count), min3, max3, (sum4 / count), min4, max4
                    }} else {{
                        print "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"
                    }}
                }}' "$CSV_FILE"
            )

            # Print to run log
            echo "Processing file: $CSV_FILE"
            echo "Mean_Weir_and_Cockerham_Fst: $Mean_Weir_and_Cockerham_Fst"
            echo "Min_Weir_and_Cockerham_Fst: $Min_Weir_and_Cockerham_Fst"
            echo "Max_Weir_and_Cockerham_Fst: $Max_Weir_and_Cockerham_Fst"
            echo "Mean_Mean_Negatives_Set_to_Zero: $Mean_Mean_Negatives_Set_to_Zero"
            echo "Min_Mean_Negatives_Set_to_Zero: $Min_Mean_Negatives_Set_to_Zero"
            echo "Max_Mean_Negatives_Set_to_Zero: $Max_Mean_Negatives_Set_to_Zero"
            echo "Mean_Mean_Threshold_Exclusion: $Mean_Mean_Threshold_Exclusion"
            echo "Min_Mean_Threshold_Exclusion: $Min_Mean_Threshold_Exclusion"
            echo "Max_Mean_Threshold_Exclusion: $Max_Mean_Threshold_Exclusion"

            # Print the data in the CSV file
            echo "$no_indv,$Mean_Weir_and_Cockerham_Fst,$Min_Weir_and_Cockerham_Fst,$Max_Weir_and_Cockerham_Fst,$Mean_Mean_Negatives_Set_to_Zero,$Min_Mean_Negatives_Set_to_Zero,$Max_Mean_Negatives_Set_to_Zero,$Mean_Mean_Threshold_Exclusion,$Min_Mean_Threshold_Exclusion,$Max_Mean_Threshold_Exclusion" >> "$Means_FILE"
        done 

        echo "Means calculation completed. Check the output in the output files within {means_path}."
    """.format(
        Output_path=Output_path,
        no_indv_list=" ".join(map(str, no_indv_list)),  # Format no_indv_list as space-separated values
        Populations=Populations,
        Input_path=Input_path,
        means_path=means_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def Get_FST_matrix_and_visualize(FSTs: list, means_path: str, CSV_path: str):
    """Extract vcftools weighted mean, max and min of each replication run and calculate the mean we are interested in."""

    inputs = {
        'FST_logs': FSTs
    }
    outputs = {
        'mean_Fst_csv': f"{means_path}/concatenated_means_all_populations.csv",
        'FST_csv': f"{CSV_path}/concatenated_runs_all_populations_all_replications.csv"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    
    # Define the shell script to be executed
    spec = """
    # Define the directory where the CSV files are located
    INPUT_DIR={means_path}
    OUTPUT_FILE={means_path}/concatenated_means_all_populations.csv

    # Initialize a flag to check if it's the first file
    first_file=true

    # Loop through each CSV file in the directory
    for file in "$INPUT_DIR"/Mean_*.csv; do
        # Extract the population identifier from the filename
        # Example: Mean_fst_SPA_TEN.csv -> SPA_TEN
        population=$(basename "$file" | awk -F'_' '{{print $3 "_" $4}}' | sed 's/.csv//')

        if [ "$first_file" = true ]; then
            # Print header for the output file using the header from the first file
            echo "population,$(head -n 1 "$file")" > "$OUTPUT_FILE"
            first_file=false  # Set the flag to false after the first file is processed
        fi

        # Read the CSV file and append data to the output file
        # Skip the header (first line) of each file, and add the population identifier as the first column
        awk -v pop="$population" 'NR > 1 {{print pop "," $0}}' "$file" >> "$OUTPUT_FILE"
    done
    echo "Concatenation completed. Output saved in $OUTPUT_FILE."

    #############NEXT CONCAT############

    # Define the directory where the CSV files are located 
    INPUT_DIR={CSV_path}
    OUTPUT_FILE={CSV_path}/concatenated_runs_all_populations_all_replications.csv

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
        means_path=means_path,
        CSV_path=CSV_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





# 2rd Job: Compute FST for files
def Calc_FST_hypo_files_plink(no_miss_hypo_vcf: str, Indvi_txt: str, Temp_path: str, Output_path: str, script_path: str, no_indv: int, replicate: int, REP_INFO: str):
    """Compute FST for subsampled individuals from the population + calculate means + Print logs"""
    inputs = {
        'no_miss_hypo_vcf': no_miss_hypo_vcf,
        'Indvi_txt': Indvi_txt
    }
    outputs = {
        'FST_file': f"{Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.FST"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    spec = """
    #Read in 'Indvi_txt' to split it without writing files.
    # Get the total number of lines in your input file
    input_file={Indvi_txt}
    total_lines=$(wc -l < "$input_file")
    half=$((total_lines / 2))

    # Create a temporary file to store the output with population labels
    output_file="{Temp_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}populations.txt"

    # Process the file to add population labels
    awk -v half="$half" 'NR <= half {{print "FAM001", $0, "population1"}} NR > half {{print "FAM002", $0, "population2"}}' "$input_file" > "$output_file"

    # Run the plink software for Fst
    plink --vcf {no_miss_hypo_vcf} --fst --within {Temp_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}populations.txt --allow-extra-chr  --double-id --set-missing-var-ids @:#\$1,\$2 --out {Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}
    #Not used yet
    #FST={Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.FST

    """.format(
        Temp_path=Temp_path, 
        Output_path=Output_path, 
        no_indv=no_indv, 
        REP_INFO=REP_INFO, 
        replicate=replicate,
        script_path=script_path,
        Indvi_txt=Indvi_txt, 
        no_miss_hypo_vcf=no_miss_hypo_vcf
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)























































# def Get_FST_means(Output_path: str, FSTs: list, no_indv_list: list, Populations: list, Input_path: str, means_path: str):
#     """Extract vcftools weighted mean, max and min of each replication run and calculate the mean we are interested in."""

#     inputs = {
#         'FST_logs': FSTs
#     }
#     outputs = {
#         'Fst_csv': f"{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_{no_indv_list[0]}.csv",
#         'Means_file':f"{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
#     }
#     options = {
#         'cores': 1,
#         'memory': '24g',
#         'walltime': '02:00:00'
#     }
    
#     # Define the shell script to be executed
#     spec = """
#         # Set directory containing FST files
#         FST_DIR="{Input_path}"
        
#         # Loop over each number in the no_indv_list
#         for no_indv in {no_indv_list}; do
#             # Define output file for this specific number of individuals
#             OUTPUT_FILE="{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_$no_indv.csv"

#             # Print the header for the CSV file
#             echo "Replicate_no,Weir_and_Cockerham_Fst,Mean_Negatives_Set_to_Zero,Mean_Threshold_Exclusion" > "$OUTPUT_FILE"
            
#             # Loop over all .fst log files matching the pattern
#             for FILE in "$FST_DIR"/*_{Populations[0]}_{Populations[1]}_FST_"$no_indv".weir.fstFST.log; do
#                 # Extract replicate number from filename (assuming pattern 'ReP_<Replicate_no>')
#                 Replicate_no=$(basename "$FILE" | grep -oP 'ReP_\\K\\d+')
                
#                 # Extract FST statistics from the file
#                 Weir_Fst=$(grep "Weir and Cockerham weighted Fst estimate" "$FILE" | awk -F': ' '{{print $2}}')
#                 Mean_Set_to_Zero=$(grep "Mean - All negative valoues are set to zero" "$FILE" | awk -F': ' '{{print $2}}')
#                 Mean_Threshold_Exclusion=$(grep "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean" "$FILE" | awk -F': ' '{{print $2}}')
                
#                 # Append the extracted data to the output CSV
#                 echo "$Replicate_no,$Weir_Fst,$Mean_Set_to_Zero,$Mean_Threshold_Exclusion" >> "$OUTPUT_FILE"
#             done

#             # Sort the CSV file by Replicate_no (first column), excluding the header
#             awk 'BEGIN{{FS = OFS = ","}} {{if (NR==1) {{print $0; next}} print $0 | "sort -t, -k1,1n "}}' "$OUTPUT_FILE" > "$OUTPUT_FILE.sorted" && mv "$OUTPUT_FILE.sorted" "$OUTPUT_FILE"
#         done  
#         echo "Data extraction completed. Check the output in the output files within {Output_path}."

#         # Print the data in the CSV file
#         echo "no_of_indv,Mean_Weir_and_Cockerham_Fst,Min_Weir_and_Cockerham_Fst,Max_Weir_and_Cockerham_Fst,Min_Mean_Negatives_Set_to_Zero,Max_Mean_Negatives_Set_to_Zero,Mean_Mean_Negatives_Set_to_Zero,Mean_Mean_Threshold_Exclusion,Min_Mean_Threshold_Exclusion,Max_Mean_Threshold_Exclusion" > "{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
#         Means_FILE="{means_path}/Mean_fst_{Populations[0]}_{Populations[1]}.csv"
#         # Loop over each number in the no_indv_list
#         for no_indv in {no_indv_list}; do
#             CSV_FILE="{Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_$no_indv.csv"
#             #Make sure premissions are given
#             chmod 777 "$CSV_FILE"
#             echo "premission given" 
#             #Calculation the mean, min and max  of the replicates - Store in variables 
#             read Mean_Weir_and_Cockerham_Fst Mean_Mean_Negatives_Set_to_Zero Mean_Mean_Threshold_Exclusion < <(awk -F',' '{{
#                 sum2 += $2; sum3 += $3; sum4 += $4; count += 1
#             }} END {{
#                 if (count > 0) {{
#                     print (sum2 / count), (sum3 / count), (sum4 / count)
#                 }} else {{
#                     print "NaN", "NaN", "NaN"
#                 }}
#             }}' "$CSV_FILE")
#             #Print to run log

#             echo "Processing file: $CSV_FILE"
#             echo "Mean_Weir_and_Cockerham_Fst: $Mean_Weir_and_Cockerham_Fst"
#             echo "Mean_Mean_Negatives_Set_to_Zero: $Mean_Mean_Negatives_Set_to_Zero"
#             echo "Mean_Mean_Threshold_Exclusion: $Mean_Mean_Threshold_Exclusion"
#             # Print the data in the CSV file
#             echo "$no_indv,$Mean_Weir_and_Cockerham_Fst,$Mean_Mean_Negatives_Set_to_Zero,$Mean_Mean_Threshold_Exclusion" >> "$Means_FILE"
#         done 

#         echo "Means calculation completed. Check the output in the output files within {means_path}."
#     """.format(
#         Output_path=Output_path,
#         no_indv_list=" ".join(map(str, no_indv_list)),  # Format no_indv_list as space-separated values
#         Populations=Populations,
#         Input_path=Input_path,
#         means_path=means_path
#     )
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)












#############A Version of  



# # 4rd Job: Extract vcftools weigthed mean, and calculate the mean we are interested in.
# def Get_FST_means(Output_path: str,  FSTs: list,  no_indv_list: list, Populations: list, Input_path: str ):
#     """Extract vcftools weigthed mean, and calculate the mean we are interested in."""
    
#     inputs = {
#         'FST_logs': FSTs
#     }
#     outputs = {
#         'Fst_csv': f"{Output_path}/complied_{Populations[0]}_{Populations[1]}.csv"
#     }
#     options = {
#         'cores': 1,
#         'memory': '24g',
#         'walltime': '02:00:00'
#     }
#     spec = """
#         #Define the directory containing your .fst files
#         FST_DIR={Input_path}

#         # Print the header to the output file
#         echo "Replicate_no,Weir_and_Cockerham_Fst,Mean_Negatives_Set_to_Zero,Mean_Threshold_Exclusion" > "$OUTPUT_FILE"

#         #Loop over the numbers in no_indiv_list
#         for no_indv in {no_indv_list}; do
#             # Define the output file
#             OUTPUT_FILE={Output_path}/complied_{Populations[0]}_{Populations[1]}_indiv_"$no_indv".csv

#             # Loop over all .fst log files in the specified directory
#             for FILE in "$FST_DIR"/*_{Populations[0]}_{Populations[1]}_FST_$no_indv".weir.fstFST.log; do
#                 # Extract the replicate number from the filename (assuming the filename has 'ReP_<Replicate_no>' pattern)
#                 Replicate_no=$(basename "$FILE" | grep -oP 'ReP_\K\d+')

#                 # Extract the required values from the file
#                 Weir_Fst=$(grep "Weir and Cockerham weighted Fst estimate" "$FILE" | awk -F': ' '{print $2}')
#                 Mean_Set_to_Zero=$(grep "Mean - All negative valoues are set to zero" "$FILE" | awk -F': ' '{print $2}')
#                 Mean_Threshold_Exclusion=$(grep "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean" "$FILE" | awk -F': ' '{print $2}')

#                 # Append the extracted values to the output file
#                 echo "$Replicate_no,$Weir_Fst,$Mean_Set_to_Zero,$Mean_Threshold_Exclusion" >> "$OUTPUT_FILE"
#             done

#         echo "Data extraction completed. Check the output in $OUTPUT_FILE."
#     """.format(
#             Output_path=Output_path,
#             no_indv_list=no_indv_list,
#             Populations=Populations, 
#             Input_path=Input_path
#     )
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





# # 1st Job: Randomly draw individuals for each population
# def Random_draw_indv(VCF: str, Temp_path: str, Populations: list, no_indv: int, Txt_path:str,):
#     """Draw a defined number of individuals from VCF from two populations."""
#     inputs = {'genome_vcf': VCF}
#     outputs = {
#         'Temp_vcf_1': f"{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz",
#         'Temp_vcf_0': f"{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz",
#         'Indvi_txt_1': f"{Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt",
#         'Indvi_txt_0': f"{Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt"
#     }
#     options = {
#         'cores': 1,
#         'memory': '16g',
#         'walltime': '01:00:00'
#     }
#     spec = f"""
#     mkdir -p {Temp_path} #Create dir, if exsist - skip

#     for pop in {' '.join(Populations)}; do

#         vcftools --vcf {VCF} --keep {Txt_path}/pop"${pop}".txt \\
#         --max-indv {no_indv} \\
#         --out {Temp_path}/"${pop}"_Tempory_{no_indv} --recode
        
#         bcftools sort {Temp_path}/"${pop}"_Tempory_{no_indv}.recode.vcf -Ob -o {Temp_path}/"${pop}"_Tempory_{no_indv}.recode.bcf.gz
        
#         bcftools index {Temp_path}/"${pop}"_Tempory_{no_indv}.recode.bcf.gz

#         bcftools query -l {Temp_path}/"${pop}"_Tempory_{no_indv}.recode.bcf.gz \\
#         > {Temp_path}/"${pop}"_Tempory_{no_indv}.txt
#     done

#     # Rename files to include both population names
#     mv {Temp_path}/{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz
#     mv {Temp_path}/{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz
#     mv {Temp_path}/{Populations[1]}_Tempory_{no_indv}.txt {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.txt
#     mv {Temp_path}/{Populations[0]}_Tempory_{no_indv}.txt {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.txt
#     mv {Temp_path}/{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz.csi {Temp_path}/{Populations[1]}_{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz.csi
#     mv {Temp_path}/{Populations[0]}_Tempory_{no_indv}.recode.bcf.gz.csi {Temp_path}/{Populations[0]}_{Populations[1]}_Tempory_{no_indv}.recode.bcf.gz.csi
#     """#.format(VCF=VCF, Temp_path=Temp_path, Populations=' '.join(Populations), no_indv=no_indv)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)











# # #Templates for RUNNING Fst on a random sampled no of indivuduals. 


# # #First job

# # # For each population, randomly draw between no_indv_list individuals!!
# # def Random_draw_indv(VCF: str, Temp_path: str, Populations: list, no_indv: int ):
# #     """Draw a defined no. of individuals from VCF From two populations"""
# #     inputs={'genome_vcf': VCF}
# #     outputs={'Temp_vcf_1': ({Temp_path}{pop}_Tempory_{no_indv}.recode.vcf.gz).format(Temp_path=Temp_path, no_indv=no_indv, Pop=${Populations[1]},
# #              'Temp_vcf_0': ({Temp_path}{pop}_Tempory_{no_indv}.recode.vcf.gz).format(Temp_path=Temp_path, no_indv=no_indv, Pop=${Populations[0]},
# #              'Indvi_txt_1':({Temp_path}{pop}_Tempory_{no_indv}.txt).format(Temp_path=Temp_path, no_indv=no_indv, Pop=${Populations[1]},
# #              'Indvi_txt_0':({Temp_path}{pop}_Tempory_{no_indv}.txt).format(Temp_path=Temp_path, no_indv=no_indv, Pop=${Populations[0]},
# #     }
# #     options={
# #     'cores': 1,
# #     'memory': '8g',
# #     'walltime': '24:00:00'
# #     }
# #     spec="""
# #         # Loop through populations
# #     for pop in "${Populations[@]}"; do
# #         # Subsample individuals and create a temporary VCF for each population
# #         vcftools --vcf {VCF} --keep ${Txt_path}"pop"${pop}.txt \
# #         --max-indv ${no_indv} \
# #         --out ${Temp_path}${pop}_Tempory_${no_indv} --recode 
# #         bgzip ${Temp_path}${pop}_Tempory_${no_indv}.vcf

# #         # Index the subsampled VCF file
# #         bcftools index ${Temp_path}${pop}_Tempory_${no_indv}.recode.vcf.gz
# #         bcftools query -l ${Temp_path}${pop}_Tempory_${no_indv}.recode.vcf.gz > ${Temp_path}${pop}_Tempory_${no_indv}.txt
# #     done
# #     """.format(VCF=vcf, Temp_path=Temp_path, Populations=','.join(Populations))
# #     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# # # 2 job  

# # # Loop over the subsampled individual (no_indv_list) and merge files !!
# # def Merge_indv_files(VCF: str, Dep_1=Dep_1, Dep_2=Dep_2, Temp_path: str, Populations: list, no_indv: int ):
# #     """ Loop over the subsampled individual files_(no_indv_list) and merge files"""
# #     inputs={'Indvi_txt_1': Dep_1,
# #             'Indvi_txt_0'Dep_2}
# #     outputs={'Merged_vcf': {Temp_path}"${Populations[0]}_${Populations[1]}"_Merged_${no_indv}.vcf.gz,
# #              'Index_csi':({Temp_path}{pop}_Tempory_{no_indv}.txt).format(Temp_path=Temp_path, no_indv=no_indv, Pop=${Populations[1]},
# #              }
# #     options={
# #     'cores': 1,
# #     'memory': '8g',
# #     'walltime': '24:00:00'
# #     }
# #     spec="""
# #         # Loop over the subsampled individual files 
# #     for no_indv in "${no_indv_list[@]}"; do

# #         # Merge the subsampled VCFs from both populations using bcftools
# #         bcftools merge \
# #         ${Temp_path}"${Populations[0]}"_Tempory_${no_indv}.recode.vcf.gz \
# #         ${Temp_path}"${Populations[1]}"_Tempory_${no_indv}.recode.vcf.gz \
# #         --threads 12 \
# #         --output-type z \
# #         --output ${Temp_path}"${Populations[0]}_${Populations[1]}"_Merged_${no_indv}.vcf.gz

# #         #Index the merged VCF file 
# #         bcftools index\
# #         --threads 6\
# #         ${Temp_path}"${Populations[0]}_${Populations[1]}"_Merged_${no_indv}.vcf.gz
# #     done
# #     """.format(VCF=vcf, Temp_path=Temp_path, Populations=','.join(Popolations))
# #     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# # # 3 job 
# # # # Loop over the subsampled individuals paird population files to compute FST
# # def Calc_FST_files(VCF: str, Temp_path: str, Output_path: str, Populations: list, no_indv: int, Replicate_no int, REP_INFO:str):
# #     """ Loop over the subsampled individuals paird population files to compute FST"""
# #     inputs={'Merged_vcf': {Temp_path}"${Populations[0]}_${Populations[1]}"_Merged_${no_indv}.vcf.gz},
# #             'Index_csi':({Temp_path}{pop}_Tempory_{no_indv}.txt).format(Temp_path=Temp_path, no_indv=no_indv, Pop=(last element of list ),
# #     outputs={'FST_file': ${Out_path}"${Populations[0]}_${Populations[1]}_FST_${no_indv}",
# #              'Index_csi':({Temp_path}{pop}_Tempory_{no_indv}.txt).format(Temp_path=Temp_path, no_indv=no_indv, Pop=(last element of list ),
# #              }
# #     options={
# #     'cores': 1,
# #     'memory': '8g',
# #     'walltime': '24:00:00'
# #     }
# #     spec="""

# #     for no_indv in "${no_indv_list[@]}"; do

# #         vcftools --gzvcf ${Temp_path}"${Populations[0]}_${Populations[1]}"_Merged_${no_indv}.vcf.gz \
# #         --weir-fst-pop ${Temp_path}"${Populations[0]}"_Tempory_${no_indv}.txt \
# #         ${Temp_path}"${Populations[1]}"_Tempory_${no_indv}.txt \
# #         --out ${Out_path}"${Populations[0]}_${Populations[1]}_FST_${no_indv}"

# #         # Check for FST computation success
# #         if [ $? -eq 0 ]; then
# #             echo  "Replecate {REP_INFO} Replicate number {REP_NO}FST calculation completed for ${no_indv} individuals in ${Populations[0]}_${Populations[1]}." >> FST.log 
# #             echo "FST calculation completed for ${no_indv} individuals in ${Populations[0]}_${Populations[1]}." >> FST.log 
# #             #echo "number of positions not nan =" (Ascript to print the #) FST.log
# #         else
# #             echo "ERROR IN CALCULAITION for ${no_indv} individuals in ${Populations[0]}_${Populations[1]}." >> FST.log 
# #         fi
# #     done

# #     """.format(VCF=vcf,Output_path=Output_path, Temp_path=Temp_path, Populations=','.join(Populations), REP_INFO=REP_INFO, REP_NO=Replicate_no)
# #     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
