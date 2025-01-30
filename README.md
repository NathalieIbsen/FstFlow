# FstFlow
Estimating and analyzing Fst values with a smooth, flexible workflow. 




Manual 1: Exploring the Fst Space cross samples
Overview
This workflow facilitates the exploration of genetic differentiation (Fst) between randomly sampled individuals across all samples in a VCF. It allows flexibility in specifying sample-size(s) and replication level. It provides detailed summaries and logs of Fst statistics. The output is a .CSV file with all replicates complied. The number of sites and number of non-variant sites are recorded in the compiled .csv file and in each run log. Three measures of Fst are outputted: 
1) The Wier Cockerham weighted mean Fst – directly from Vcftools. 
2) A mean Fst where all negative values are excluded. 
3) A mean Fst where all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded. 
The steps of the workflow are detailed at the end of the manual.   
Input Requirements

A filtered .vcf file containing all samples.
The .vcf can include either all sites, or only SNPs. Including all sites does not affect Fst estimation, as monomorphic sites will result in “-nan”, but may heavily impact computational speed, especially is using whole-genome data.
Ensure the filtering matches the type of sequencing data and current recommendations.
Running the Workflow
Prepare the Environment
Install the required software and compatible versions using the ENV.yaml file.
e.g. through conda 
(script_folder)$ conda create -n EVN  
Prepare the workflow.py:
Place the workflow.py and templates.py in a folder. From here the pipeline will run, unless otherwise manually specified.
Declare the required input files and workflow variables:
#########################
# Replication information 
#########################
# Replication information (metadata reference)
REP_INFO = 'First_run'
# Number of FULL RUN replications to run
no_replicates = 1000
# List of how many individuals to sample for each hypothetical population comparison
no_indv_list = [4, 6, 8, 10]
#OBS! Each value corresponds to a number of individuals that be sampled before they are spilt into populations. Eg if 4 then 2 individuals will be sampled from each hypothetical population. If 10 then 5 individuals will be sampled from each population and so on. 


# can also be a single int, input as [int]  


##############
# INPUT FILES
##############
#Population vcf 
VCF='/Full/Path/To/.vcf'



Specify the path to the script folder containing the workflow.py & templates.py:
##############################################
# PATHS FOR workflow OUTPUTs, RESULTS AND TEMP FILES
##############################################
#Where the workflow and templates file are:
scripts_path='/Full/Path/To/scripts_folder'



Optionally:
Input paths manually.  
Optionally: 
Choose to not “clean” the workflow. This will keep all temporary files.
##############
#OPTIONAL  !!!
##############
# If you want the temporary files, VCFs, pr site FST estimates and the log files to remain after your run input NO in the below variable. 


clean_workflow='YES' #optinally 'YES' or 'NO' 



Execute the Workflow
Navigate to the script folder and execute the following commands:
(script_folder)$ gwf run          # Send the jobs to the queue  
(script_folder)$ gwf status       # Display the status of each job  
(script_folder)$ gwf status -f summary #Will display a count of each status within workflow
(script_folder)$ gwf logs "$Job_ID"  # Print std.out logs for a specific job
(script_folder)$ gwf logs -e "$Job_ID"  # Print ERROR logs for a specific job
(script_folder)$ gwf info "$Job_ID"  #will print the bash script as interpreted by gwf.
Please see the gwf documentation for detailed information on usage and even cooler utilities.

The output of the pipeline
By default:
In the RESULTS folder three folders will be created called 1) “$no_replicates”_replicates, 2) CSV, 3) MEAN. 
Within 1) “$no_replicates”_replicates folder the *.wier.fst output from Vcftools and the created  *FST.log will be placed. By default, the *.wier.fst output is deleted upon completion of the workflow. - If you wish to keep these, please change the clean_workflow parameter to “NO”, in the workflow. 
Within 2) CSV folder the complied Fst estimates and metrics for each run at a set sample-size. There will be a .csv file for each sample-size used. The .csv contains the following headers and data: 
Replicate_no, = a numeric value
Weir_and_Cockerham_Fst, = The Wier Cockerham weighted mean Fst – directly from Vcftools.
Mean_Negatives_Set_to_Zero, = A mean Fst where all negative values are excluded.
Mean_Threshold_Exclusion, = A mean Fst where all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded.
no_sites , = number of sites in which genotypes were called for all samples. 
no_monomorph = number of monomorphic sites cross the samples. 
The number of SNPs analyzed can be calculated = no_sites - no_monomorph
Within 3) MEAN folder a .csv file containing means and metrics cross a sample-size for all sample-sizes used. I will contain the following headers and data: 
no_of_indv, = number of individuals drawn
Mean_Weir_and_Cockerham_Fst, = The mean weighted Wier Cockerham Fst 
Min_Weir_and_Cockerham_Fst, = The minimum mean weigthed Wier Cockerham Fst
Max_Weir_and_Cockerham_Fst, = The Maximum mean weigthed Wier Cockerham Fst         
Mean_Mean_Negatives_Set_to_Zero, = The mean of all means when Negatives are set to zero
Min_Mean_Negatives_Set_to_Zero, = The lowest mean when Negatives are set to zero
Max_Mean_Negatives_Set_to_Zero, = The highest mean when Negatives are set to zero              
Mean_Mean_Threshold_Exclusion, = The mean of all means when all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded.
Min_Mean_Threshold_Exclusion, = The lowest mean when all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded.
Max_Mean_Threshold_Exclusion, = The highest mean when all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded.

Steps in the Workflow
Subset Creation
Generate subsets of individuals from a filtered .vcf file to simulate hypothetical populations.
Randomly drawing individuals to simulate hypothetical populations (Vcftools V0.1.16)
Merged, sort, index and extract sample IDs (Bcftools V1.19)
The .vcf file is filtered for missing data, leaving only positions with genotypes called for all samples. (Vcftools V0.1.16)
Fst Computation
Compute the Fst between two subsets of individuals within the hypothetical population.
This step is repeated for the number of replications as specified. (Vcftools V0.1.16)
Post each estimation a custom bash script (standard bash) is used to calculate metrics of the run and extract Fst values. 
Result Compilation and Summarization
Compile and summarize Fst results across all replicates to produce the summary of genetic differentiation. (Custom script – standard bash). 




Manual 2: Estimating and Summarizing Fst
Overview
This workflow generates all unique combinations of individuals for each pairwise population comparison at specified sample sizes. It allows flexibility in specifying sample-size(s) and replication level. It provides detailed summaries and logs of Fst statistics. The output is a .CSV file with all replicates complied. The number of sites and number of non-variant sites are recorded in the compiled .csv file and in each run log. Three measures of Fst are outputted: 
1) The Wier Cockerham weighted mean Fst – directly from Vcftools. 
2) A mean Fst where all negative values are excluded. 
3) A mean Fst where all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded. The steps of the workflow are detailed at the end of the manual.   
Input Requirements
A filtered .vcf file containing all samples.
The .vcf can include either all sites, or only SNPs. Including all sites does not affect Fst estimation, as monomorphic sites will result in “-nan”, but may heavily impact computational speed, especially is using whole-genome data.
Ensure the filtering matches the type of sequencing data and current recommendations.
A set of population files (pop***.txt) for each population.
***= a 3 digit population identifier. Files must contain one sample name each line and adhere to the naming of samples in the .VCF file.
Running the Workflow
Prepare the Environment
Install the required software and compatible versions using the ENV.yaml file.
e.g. through conda 
(script_folder)$ conda create -n EVN  
Prepare the workflow.py:
Place the workflow.py and templates.py in a folder. From here the pipeline will run, unless otherwise manually specified.
Declare the required input files: 
##############
# INPUT FILES#
##############


#Population vcf 
VCF='/Full/Path/To/.vcf'


#Path to a folder where pop***.txt files are located (***=3 digit population reference)
Txt_path = '/full/path/to/Population_lists'



 and workflow variables:
#########################################
## Replication INFORMATION and OPTIONS ##
#########################################


# Replication information 
REP_INFO = 'min_samples_2'  # Replication information (metadata reference .log REFerence)
#Minimun number of samples pr. population to be consindered
min_samples = 2  #Input (numeric: int) 
# List of how many individuals to sample from each population 
no_indv_list = [2]
#OBS! ^Maxmium number in no_indv_list must be =<min_samples. [int] or [L,I,S,T] but always protected by [] 



Specify the path to the script folder containing the workflow.py & templates.py:
##############################################
# PATHS FOR workflow OUTPUTs, RESULTS AND TEMP FILES
##############################################
#Where the workflow and templates file are:
scripts_path='/Full/Path/To/scripts_folder'



Optionally:
Input paths manually.  
Optionally: 
Choose to not “clean” the workflow. This will keep all temporary files.
##############
#OPTIONAL  !!!
##############
# If you want the temporary files, VCFs, pr site FST estimates and the log files to remain after your run input NO in the below variable. 


clean_workflow='YES' #optinally 'YES' or 'NO' 




Execute the Workflow
Navigate to the script folder and execute the following commands:
(script_folder)$ gwf run          # Send the jobs to the queue  
(script_folder)$ gwf status       # Display the status of each job  
(script_folder)$ gwf status -f summary #Will display count of each status with in workflow
(script_folder)$ gwf logs "$Job_ID"  # Print std.out logs for a specific job
(script_folder)$ gwf logs -e "$Job_ID"  # Print ERROR logs for a specific job
(script_folder)$ gwf info "$Job_ID"  #will print the bash script as interpreted by gwf.
Please see the gwf documentation for detailed information on usage and even cooler utilities.

The output of the pipeline
By default:
In the RESULTS folder two folders will be created called 1) FST, 2) CSV. 
Within 1) FST folder the *.wier.fst output from Vcftools and the created *FST.log will be placed. By default, this folder is deleted upon completion of the workflow. - If you wish to keep these files, please change the clean_workflow parameter to “NO”, in the workflow. 
Within 2) CSV folder the complied Fst estimates and metrics for each run at a set sample size. There will be a .csv file for each population pair and sample size used, as well as a compiled .csv. These files contain the following headers and data: 
Replicate_no, = a numeric value
Weir_and_Cockerham_Fst, = The Wier Cockerham weighted mean Fst – directly from Vcftools.
Mean_Negatives_Set_to_Zero, = A mean Fst where all negative values are excluded.
Mean_Threshold_Exclusion, = A mean Fst where all Fst’s below 0 but above -0.10 is treated as a 0, and Fst values below -0.1 are excluded.
no_sites , = number of sites in which genotypes were called for all samples. 
no_monomorph = number of monomorphic sites across the samples. 
The complied .csv will also contain a collum called no_of_indv, = number of individuals drawn.
The number of SNPs analyzed can be calculated = no_sites - no_monomorph

Steps in the Workflow
Generate Unique Combinations
Create all unique combinations of individuals for each population pair based on the provided .vcf and pop***.txt files.
Population files should have a three-digit population identifier and follow the naming convention pop***.txt.
Fst Estimation
Process the .vcf file and estimate Fst values for each unique combination of population pairs.
Perform calculations for specified sample sizes and multiple replicates.
Result Compilation
Compile results from each population pair and replicate.
Consolidate all runs into a single .CSV file for ease of downstream analysis.


