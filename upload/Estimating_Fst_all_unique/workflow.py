from gwf import Workflow, AnonymousTarget
import os
import itertools
from itertools import combinations
from templates import *  # Assuming you have templates defined elsewhere

gwf = Workflow(defaults={'account': 'EcoGenetics'})

# ________________________________________________________________________________________
#                             ##############################
#                             ## Input variables for the RUN 
#                             ##############################

##############
# INPUT FILES#
##############

#Population vcf 
VCF = '/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/scripts/sorted_Cyrtophora_full4.vcf'
#Path to a folder where pop***.txt files are located (***=3 digit population reference)
Txt_path = '/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/Population_lists'

#########################################
## Replication INFORMATION and OPTIONS ##
#########################################

# Replication information 
REP_INFO = 'min_samples_2'  # Replication information (metadata reference .log REFerence)
#Minimun number of samples pr. population to be consindered
min_samples = 2  #Input (numeric: int) 
# List of how many individuals to sample from each population 
no_indv_list = [2]
#OBS! ^Maxmium number in no_indv_list must be =<min_samples. int or list  

####################################################
# PATHS FOR workflow OUTPUTs, RESULTS AND TEMP FILES
####################################################

#Where the workflow and templatesfile is:
scripts_path='/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/scripts'


#_________________________________________OPTIONS_____________________________________________

##############
#OPTIONAL  !!!
##############
# If you want the temporary files, VCFs, pr site FST estimates and the log files to remain after your run input NO in the below variable. 

clean_workflow='YES' #optinally 'YES' or 'NO' 

##############
#OPTIONAL !!! 
##############
# Manually input the paths to where you would like the files to be located. 
#If these folders do not exsist, the workflow will create them.

#Raw vcftools pr. site estimates 
Out_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/RESULTS/All_unique/FST"
#Temp folder 
Temp_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/temps/All_unique"
#Results Folder
CSV_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/RESULTS/All_unique/CSV"


# #Raw vcftools pr. site estimates 
# Out_path = f"{scripts_path}/{REP_INFO}/RESULTS/{REP_INFO}/{no_replicates}_replicates"
# #Temp folder 
# Temp_path = f"{scripts_path}/{REP_INFO}/temps/{REP_INFO}/{no_replicates}_replicates"
# #Results Folder
# CSV_path = f"{scripts_path}/{REP_INFO}/RESULTS/{REP_INFO}/CSV"



##########################################
# Create the folders if they don't exist #

os.makedirs(Out_path, exist_ok=True)
os.makedirs(Temp_path, exist_ok=True)
os.makedirs(CSV_path, exist_ok=True)

print(f"Created or confirmed: {Out_path}")
print(f"Created or confirmed: {Temp_path}")
print(f"Created or confirmed: {CSV_path}")

#________________________________________________________________________________#
#DON'T CHANGE!!

##############C########################
# Variables created by the work-flow #
######################################

Populations_list = []  # List to store population identifiers

####ONLY populations with at least 5 inviduals >=min_sa
# Iterate over files in the population lists directory
for filename in os.listdir(Txt_path):
    file_path = os.path.join(Txt_path, filename)
    # Ensure it's a file before proceeding
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            # Only add populations with at least 'min_samples' individuals
            if len(lines) >= min_samples:
                population_id = filename[3:6]  # Extract characters 4-6 from the filename
                Populations_list.append(population_id)


# Generate all pairwise combinations of population IDs
Populations = list(combinations(Populations_list, 2))
print("Population list done")

#But... 
#debug
#IF you want it printed: un-#
print("Populations List:", Populations_list)


#If You want it printed:  un-
# print("Pairwise Population IDs:")
# for pair in Populations:
#     print(pair)


# ________________________________________________________________________________________
#                             ###########################
#                                  Jobs in Work-flow  
#                             ###########################


all_FST_CSV = []

#"""""""""""""""""""""""""""""""""""
###FOR later input option !!!!
#""""""""""""""""""""""""""""""""""
# Iterate over individual counts and population pairs to define workflow jobs

for i, pop in enumerate(Populations):
        #Draw a defined number of individuals from VCF from two populations.
        # + sort & zip vcf 
        # + Index 
        # + Log individual files (metadata): ->  "{Temp_path}/${{pop}}_Tempory_{no_indv}.txt"
            job_name = f"RDI_" + ''.join(pop) # Define a unique job name for each individual count and population pair
            job1 = gwf.target_from_template(
                name=job_name + REP_INFO,
                template=Create_all_unique(
                    VCF=VCF,
                    Temp_path = Temp_path + '/' + ''.join(pop) + '/',
                    Populations=pop,  # Pass the population pair as a tuple
                    no_indv_list=no_indv_list,
                    Txt_path=Txt_path  + '/',
                    results_path=Out_path + '/',
                    REP_INFO=REP_INFO,
                    CSV_path=CSV_path + '/'
                )
            )
            Fst_CSV=job1.outputs['no_indv']
            all_FST_CSV.append(Fst_CSV)

#debug
#Unhash the below and print all CSV to the gwf.log
#print(', '.join(all_FST_CSV))

#compile all data
job2 = gwf.target_from_template(
    name='Get_FST_compiled_CSV',
    template=Get_FST_compiled_CSV(
        CSVs=all_FST_CSV,
        CSV_path=CSV_path
                    )
                )

#clean workflow 
if clean_workflow == "YES":
    job3 = gwf.target_from_template(
    name='Cleaning',
    template=Clean(
        CSV=job2.outputs['FST_csv'],
        Out_path=Out_path,
        CSV_path=CSV_path + '/',
        Temp_path = Temp_path
                    )
                )
else:  # If clean_workflow is not "YES"
    print("Temporary files, VCFs, pr. site Fst, and log files not removed")






































#             #     ## THIS JOB HAS BEEN ADDed TO the job one function script. 
#             # # #merge pop files     
#             # job2 = gwf.target_from_template(
#             #     name=f"Merge_Pop_{replicate}{pop[0]}{pop[1]}{no_indv}" + REP_INFO, #-"-"
#             #     template=Merge_indv_files(
#             #         Temp_path = Temp_path + '/' + ''.join(pop) + '/' + str(replicate),
#             #         Populations=pop,  # Pass the population pair as a tuple
#             #         no_indv=no_indv,
#             #         Temp_vcf_0=job1.outputs['Temp_vcf_0'],
#             #         Temp_vcf_1=job1.outputs['Temp_vcf_1']
#             #     )
#             # )
#             # #Estimate FST





























































































# #_________-old run

# # from gwf import Workflow, AnonymousTarget


# # import glob
# # from templates import *
# # from itertools import combinations 

# # gwf = Workflow(defaults={'account': 'EcoGenetics'})

# # #________________________________________________________________________________________
# #                                 ##############################
# #                                 ## Input variables for the RUN 
# #                                 ##############################


# # #REPLICATION INFORMATION 
# # #############################

# # #Replication SETTINGS !! 
# # #############

# # # Replication information (Name reference for metadata)
# # REP_INFO='First_trial'

# # ## Numnber of replications to run. (A new set of indivuduals will be sampled pr. run.)
# # no_replicates=1

# # ### List of individual counts to sample from each population
# # no_indv_list=[1, 2, 3, 4, 5]


# # #### INPUT FILES 
# # ################
# # # 
# # # Population level VCF PATH/INPUT.vcf
# # VCF='/faststorage/project/EcoGenetics/people/Tammy_Ho/RADseq/Cyrtophora_full4_outfiles/Cyrtophora_full4.vcf'

# # ##Population 
# # Txt_path='/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/Population_lists/'

# # ##########
# # ##WHERE TO STORE THINGS ?
# # ##########

# # # Output_ RESLUTS path
# # Out_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/RESULTS/{REP_INFO}/{no_replicates}_replicates"

# # #Temporary_path
# # Temp_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/temps/{REP_INFO}/{no_replicates}_replicates"


# # # Create the folders if they don't exist
# # os.makedirs(Out_path, exist_ok=True)
# # os.makedirs(Temp_path, exist_ok=True)

# # print(f"Created or confirmed: {Out_path}")
# # print(f"Created or confirmed: {Temp_path}")


# # #######Variabels created by the work-flow
# # ##########################################

# # # Initialize an empty list to store population identifiers
# # Populations_list = []

# # # Iterate through each file in the folder
# # for filename in os.listdir(Txt_path):
# #     # Ensure we only process files (ignore directories)
# #     file_path = os.path.join(Txt_path, filename)
# #     if os.path.isfile(file_path):
# #         # Count the number of lines in the file
# #         with open(file_path, 'r') as file:
# #             lines = file.readlines()
# #             if len(lines) > 4:
# #                 # Extract population identifier (characters 4-7 of the filename)
# #                 population_id = filename[3:6]
# #                 Populations_list.append(population_id)

# # print("Populations List:", Populations_list)

# # # Generate all pairwise combinations of population IDs
# # Populations = list(combinations(Populations_list, 2))

# # # Print the pairwise combinations
# # print("Pairwise Population IDs:")
# # for pair in Populations:
# #     print(pair)
# # #__________________________________________________________________________________________
# #                                     ###########################
# #                                         #Jobs in Work-flow
# #                                     ##########################

# # for i, no_indiv in enumerate(no_indv_list):
# #     for i, pop in enumerate(Populaions):
# #     # For each population, randomly draw between no_indv_list individuals!!
# #         #gwf move and rename ped and maps  
# #         job2 = gwf.target_from_template(
# #             name='Random_draw_indv' + no_replicates + pop[1] + pop[0],
# #             template=Random_draw_indv(
# #                 VCF=VCF,
# #                 Temp_path=Temp_path,
# #                 Populations=','.join(pop),
# #                 no_indv=no_indv
# #                 )
# #             )