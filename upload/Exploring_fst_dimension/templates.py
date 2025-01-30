from gwf import AnonymousTarget
import os

############################################################################
# #Templates for EXPLORING THE Fst space of all samples in a vcf#
############################################################################

#______________________________________STEP 1_____________________________________________


def Random_hypo_draw(VCF: str, Temp_path: str, replicate: str, no_indv: int):
    """1)Draw a defined number of individuals from vcf = hypothetical population.
    2)Sort, index and write a popualtion file, 
    3)Filter out missing positions,
    4) write out .txt file with population"""
    
    inputs = {'genome_vcf': VCF}
    outputs = {
        'Indvi_txt': f"{Temp_path}/{replicate}_Tempory_{no_indv}.txt",
        'no_miss_hypo_vcf': f"{Temp_path}/{replicate}no_miss_hypo_{no_indv}.bcf.gz"
    }
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '01:00:00'
    }
    spec =f"""
    mkdir -p {Temp_path}  # Create directory, if it exists - skip

    vcftools --vcf {VCF} --max-indv {no_indv} \\
    --max-missing 1 \\
    --out {Temp_path}/{replicate}_Tempory_{no_indv} --recode --recode-INFO-all

    bcftools sort {Temp_path}/{replicate}_Tempory_{no_indv}.recode.vcf -Ob \\
    -o {Temp_path}/{replicate}no_miss_hypo_{no_indv}.bcf.gz

    bcftools index {Temp_path}/{replicate}no_miss_hypo_{no_indv}.bcf.gz

    bcftools query -l {Temp_path}/{replicate}no_miss_hypo_{no_indv}.bcf.gz \\
    > {Temp_path}/{replicate}_Tempory_{no_indv}.txt
    
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


#______________________________________STEP 2_____________________________________________


# 2rd Job: Compute FST for files
def Calc_FST_vcftools(no_miss_hypo_vcf: str, Indvi_txt: str, Temp_path: str, Output_path: str, script_path: str, no_indv: int, replicate: int, REP_INFO: str):
    """Compute FST for subsampled individuals from the population + calculate means + Print logs"""
    inputs = {
        'no_miss_hypo_vcf': no_miss_hypo_vcf,
        'Indvi_txt': Indvi_txt
    }
    outputs = {
        'FST_file': f"{Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.weir.fst", 
        'log': f"{Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.weir.fstFST.log"
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

    mkdir -p {Temp_path}/Hypo_pops/
    # Create temporary files for the top and bottom halves
    first_half_file="{Temp_path}/Hypo_pops/hypo_{replicate}_{no_indv}_a.txt"
    bottom_half_file="{Temp_path}/Hypo_pops/hypo_{replicate}_{no_indv}_b.txt"

    # Extract the top half of the file
    head -n "$half" "$input_file" > "$first_half_file"

    # Extract the bottom half of the file
    tail -n "$half" "$input_file" > "$bottom_half_file"

    # Run the Fst software
    bcftools view -O v {no_miss_hypo_vcf} | vcftools --vcf - \\
        --weir-fst-pop "$first_half_file" \\
        --weir-fst-pop "$bottom_half_file" \\
        --out {Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}

    FST={Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.weir.fst

    #A fix is needed to save the .out of VCFtools  - maybe okay now ?  OBS!
    gwf_Out="{script_path}/.gwf/logs/FST_vcftools_{replicate}_{no_indv}{REP_INFO}.stderr"
    # Print out two logs# 
    if [ $? -eq 0 ]; then

        #### First log: individual FST.log#
        
        echo "Replicate {REP_INFO} - Replicate {replicate}: FST calculation completed for {no_indv} individuals." >> "$FST"FST.log
        #input sample info
        cat "$first_half_file" >> "$FST"FST.log
        cat "$bottom_half_file" >> "$FST"FST.log

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
        echo "______Mean estimates Replicate {replicate}__________" >> "$FST"FST.log
        echo "Weir and Cockerham weighted Fst estimate: $fst_value" >> "$FST"FST.log
        echo "Mean - All negative valoues are set to zero: $mean_neg_to_zero" >> "$FST"FST.log
        echo "Mean - All negative valoues >-0.10 are set to zero, any number <-0.10 are excluded from the mean: $mean_above_threshold" >> "$FST"FST.log

    else
        echo "ERROR: FST calculation failed for {no_indv} individuals." >> "$FST"FST.log
        echo "$first_half_file" >> "$FST"FST.log
        echo "$first_half_file">> "$FST"FST.log

        # Error handling: Log errors and file paths if FST calculation fails
        echo "ERROR: FST calculation failed for {no_indv} individuals." >> {Output_path}FST_RUN.log
        echo "$first_half_file" >> {Output_path}FST_RUN.log
        echo "$first_half_file" >> {Output_path}FST_RUN.log
    fi 

    # Optionally, clean up the temporary files
    rm "$first_half_file" "$bottom_half_file"
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


#______________________________________STEP 3_____________________________________________


def Get_hypo_FST_means(Output_path: str, FSTs: list, no_indv_list: list, replicate: int, Input_path: str, means_path: str):
    """Extract vcftools weighted mean, max and min of each replication run and calculate the mean we are interested in."""
    inputs = {
        'FST_logs': FSTs
    }
    outputs = {
        'Fst_csv': f"{Output_path}/complied_FST_indiv_{no_indv_list[0]}.csv",
        'Means_file':f"{means_path}/Mean_fst.csv"
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    # Define the shell script to be executed
    spec = """
        # Set directory containing FST files
        FST_DIR="{Input_path}"
        
        # Loop over each number in the no_indv_list
        for no_indv in {no_indv_list}; do
            # Define output file for this specific number of individuals
            OUTPUT_FILE="{Output_path}/complied_FST_indiv_$no_indv.csv"

            # Print the header for the CSV file
            echo "Replicate_no,Weir_and_Cockerham_Fst,Mean_Negatives_Set_to_Zero,Mean_Threshold_Exclusion,no_sites,no_monomorph" > "$OUTPUT_FILE"
            
            # Loop over all .fst log files matching the pattern
            for FILE in "$FST_DIR"/*_FST_"$no_indv".weir.fstFST.log; do
                # Extract replicate number from filename (assuming pattern 'ReP_<Replicate_no>')
                Replicate_no=$(basename "$FILE" | grep -oP 'hypo_\\K\\d+')
                
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
        echo "no_of_indv,Mean_Weir_and_Cockerham_Fst,Min_Weir_and_Cockerham_Fst,Max_Weir_and_Cockerham_Fst,\\
            Mean_Mean_Negatives_Set_to_Zero,Min_Mean_Negatives_Set_to_Zero,Max_Mean_Negatives_Set_to_Zero,\\
                Mean_Mean_Threshold_Exclusion,Min_Mean_Threshold_Exclusion,Max_Mean_Threshold_Exclusion"\\
                      > "{means_path}/Mean_fst.csv"
        Means_FILE="{means_path}/Mean_fst.csv"
        # Loop over each number in the no_indv_list
        for no_indv in {no_indv_list}; do
            CSV_FILE="{Output_path}/complied_FST_indiv_$no_indv.csv"
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
        replicate=replicate,
        Input_path=Input_path,
        means_path=means_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



#______________________________________STEP 4_____________________________________________


def Clean(CSVs: str, Temp_path: str, Out_path: str):
    """Cleaning up post run, remove all temperoary files"""
    inputs = {
        'CSVs': CSVs
    }
    outputs = {
    }
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '02:00:00'
    }
    # Define the shell script to be executed
    spec = """
    echo "Removing temporary files from {Temps_path}..."
    rm -rf {Temps_path}
    echo "Removing .fst files from {Out_path}..."
    find {Out_path} -type f -name "*.fst" -exec rm -f {{}} \;
    """.format(
        Temps_path=Temp_path,
        Out_path=Out_path
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


































































# # 2rd Job: Compute FST for files WITH PLINK!!!

# def Calc_FST_hypo_files_plink(dependency:str, no_miss_hypo_vcf: str, Indvi_txt: str, Temp_path: str, Output_path: str, script_path: str, no_indv: int, replicate: int, REP_INFO: str):
#     """Compute FST for subsampled individuals from the population + calculate means + Print logs"""
#     inputs = {
#         'no_miss_hypo_vcf': no_miss_hypo_vcf,
#         'Indvi_txt': Indvi_txt,
#         'dependency': dependency,
#     }
#     outputs = {
#         'FST_log': f"{Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.log"
#     }
#     options = {
#         'cores': 1,
#         'memory': '24g',
#         'walltime': '02:00:00'
#     }
#     spec = """
#     #Read in 'Indvi_txt' to split it without writing files.
#     # Get the total number of lines in your input file
#     input_file={Indvi_txt}
#     total_lines=$(wc -l < "$input_file")
#     half=$((total_lines / 2))

#     # Create a temporary file to store the output with population labels
#     output_file="{Temp_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}populations.txt"
#     # Process the file to add population labels
#     awk -v half="$half" 'NR <= half {{print $0, $0, "population1"}} NR > half {{print $0, $0, "population2"}}' "$input_file" > "$output_file"

#     #FIX VCF FOR PLINK
#     # input 0 for chr name and replape position with current row position.
#     bcftools view -O v {no_miss_hypo_vcf} | \
#     awk 'BEGIN {{OFS="\t"}} 
#         /^#/ {{print; next}} 
#         {{ $1 = "0"; $2 = NR - 1; print }}' > {Temp_path}/TEMP{REP_INFO}_hypo_{replicate}_FST_{no_indv}.vcf

#     echo "Filtered VCF written to {Temp_path}/TEMP{REP_INFO}_hypo_{replicate}_FST_{no_indv}.vcf"

#     # Run the plink software for Fst
#     plink --fst --vcf {Temp_path}/TEMP{REP_INFO}_hypo_{replicate}_FST_{no_indv}.vcf --within {Temp_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}populations.txt --allow-extra-chr  --double-id --set-missing-var-ids @:#\$1,\$2 --out {Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}
#     #Not used yet
#     #FST={Output_path}/{REP_INFO}_hypo_{replicate}_FST_{no_indv}.FST

#     """.format(
#         Temp_path=Temp_path, 
#         Output_path=Output_path, 
#         no_indv=no_indv, 
#         REP_INFO=REP_INFO, 
#         replicate=replicate,
#         script_path=script_path,
#         Indvi_txt=Indvi_txt, 
#         no_miss_hypo_vcf=no_miss_hypo_vcf
#     )
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)








# def Get_hypo_FST_means_plink(Output_path: str, FSTs: list, no_indv_list: list, replicate: int, Input_path: str, means_path: str):
#     """Extract plink fat mean and weighted mean."""
#     inputs = {
#         'FST_logs': FSTs
#     }
#     outputs = {
#         'Fst_csv': f"{Output_path}/plink_complied_FST_indiv_{no_indv_list[0]}.csv",
#         'Means_file_plink':f"{means_path}/Mean_fst_plink.csv"
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
#             OUTPUT_FILE="{Output_path}/plink_complied_FST_indiv_$no_indv.csv"

#             # Print the header for the CSV file
#             echo "Replicate_no,mean_Fst,weighted_mean_Fst,no_of_markers" > "$OUTPUT_FILE"
            
#             # Loop over all .fst log files matching the pattern
#             for FILE in "$FST_DIR"/*_FST_"$no_indv".log; do
#                 # Extract replicate number from filename (assuming pattern 'hypo_<Replicate_no>')
#                 Replicate_no=$(basename "$FILE" | grep -oP 'hypo_\\K\\d+')
                
#                 # Extract FST statistics from the file
#                 Weigthed_mean_Fst=$(grep "Weighted Fst estimate" "$FILE" | awk -F': ' '{{print $2}}')
#                 Mean_Fst=$(grep "Mean Fst estimate" "$FILE" | awk -F': ' '{{print $2}}') 
#                 No_of_sites=$(grep "markers with valid Fst estimates" "$FILE" | awk -F' ' '{{print $0}}')
#                 echo "$Replicate_no,$Weigthed_mean_Fst,$Mean_Fst"
#                 # Append the extracted data to the output CSV
#                 echo "$Replicate_no,$Mean_Fst,$Weigthed_mean_Fst,$No_of_sites" >> "$OUTPUT_FILE"
#             done

#             # Sort the CSV file by Replicate_no (first column), excluding the header
#             awk 'BEGIN{{FS = OFS = ","}} {{if (NR==1) {{print $0; next}} print $0 | "sort -t, -k1,1n "}}' "$OUTPUT_FILE" > "$OUTPUT_FILE.sorted" && mv "$OUTPUT_FILE.sorted" "$OUTPUT_FILE"
#         done  
#         echo "Data extraction completed. Check the output in the output files within {Output_path}."

        
#         # Print the data in the MEANS !!! CSV file
#         echo "mean_Fst,min_mean_Fst,max_mean_Fst,weighted_mean_Fst,min_weighted_mean_Fst,max_weighted_mean_Fst"\\
#                       > "{means_path}/Mean_fst_plink.csv"
#         Means_FILE="{means_path}/Mean_fst_plink.csv"
#         # Loop over each number in the no_indv_list
#         for no_indv in {no_indv_list}; do
#             CSV_FILE="{Output_path}/plink_complied_FST_indiv_$no_indv.csv"
#             #Make sure premissions are given
#             chmod 777 "$CSV_FILE"
#             # Calculate the mean, min, and max of the replicates - Store in variables 
#             read mean_Fst min_mean_Fst  max_mean_Fst weighted_mean_Fst min_weighted_mean_Fst max_weighted_mean_Fst < <(
#                 awk -F',' 'NR > 1 {{ # skip header 
#                     sum2 += $2; sum3 += $3; sum4 += $4;
#                     if ($2 < min2 || NR == 2) min2 = $2;
#                     if ($2 > max2 || NR == 2) max2 = $2;
#                     if ($3 < min3 || NR == 2) min3 = $3;
#                     if ($3 > max3 || NR == 2) max3 = $3;
#                     count += 1
#                 }} END {{
#                     if (count > 0) {{
#                         print (sum2 / count), min2, max2, (sum3 / count), min3, max3
#                     }} else {{
#                         print "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"
#                     }}
#                 }}' "$CSV_FILE"
#             )

#             # Print the data in the CSV file
#             echo "$mean_Fst,$min_mean_Fst,$max_mean_Fst,$weighted_mean_Fst,$min_weighted_mean_Fst,$max_weighted_mean_Fst" >> "$Means_FILE"
#         done 
#         echo "Means calculation completed. Check the output in the output files within {means_path}."
#     """.format(
#         Output_path=Output_path,
#         no_indv_list=" ".join(map(str, no_indv_list)),  # Format no_indv_list as space-separated values
#         replicate=replicate,
#         Input_path=Input_path,
#         means_path=means_path
#     )
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


