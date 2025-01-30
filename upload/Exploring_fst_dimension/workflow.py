from gwf import Workflow, AnonymousTarget
import os
from itertools import combinations
from templates import * 

gwf = Workflow(defaults={'account': 'EcoGenetics'})

# ________________________________________________________________________________________
#                             ##############################
#                             ## Input variables for the RUN 
#                             ##############################


#############################
## Replication INFORMATION ##
#############################

# Replication information 
REP_INFO = 'FST_without_TEN'  # Replication information (metadata reference .log REFerence)
# Number of FULL RUN replications to run
no_replicates = 1000
# List of how many individuals to sample to form new populations 4 = sample size of 2, 
# 6 = sample size of 3, and so on.  
no_indv_list = [4, 6, 8, 10]


##############
# INPUT FILES
##############

#Population vcf 
VCF = '/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/scripts/new_vcf_no_TEN.recode.vcf'

###############
# scripts path 
###############
#The path to where the templates.py and the workflow.py is located

script_path='/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/scripts/Exploring_fst_dimension'


#________________________________________OPTIONS___________________________________________

##############
#OPTIONAL!!!
##############
# If you want the temporary files, VCFs, pr site FST estimates and the log files to remain after your run input NO in the below variable. 

clean_workflow='YES' #optinally 'YES' or 'NO' 


##############
#OPTIONAL !!! 
##############
# Manually input the paths to where you would like the workfloe run and files to be located. 
#If these folders do not exsist, the workflow will create them.

Fst_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/RESULTS/{no_replicates}_replicates"
CSV_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/RESULTS/CSV"
Means_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/RESULTS/MEANS"
Temp_path = f"/faststorage/project/EcoGenetics/people/Nathalie/Cyrtophora/{REP_INFO}/temps/{no_replicates}_replicates"

# #Raw vcftools pr. site estimates 
# Out_path = f"{scripts_path}/{REP_INFO}/RESULTS/{REP_INFO}/{no_replicates}_replicates"
# #Temp folder 
# Temp_path = f"{scripts_path}/{REP_INFO}/temps/{REP_INFO}/{no_replicates}_replicates"
# #Results Folder
# CSV_path = f"{scripts_path}/{REP_INFO}/RESULTS/{REP_INFO}/CSV"
# Means folder  
# Means_path = f"{scripts_path}/{REP_INFO}/RESULTS/{REP_INFO}/MEANS"



# ___________________________________________________________________#
#DON'T CHANGE!!

######################################
# Variables created by the work-flow #

##########################################
# Create the folders if they don't exist #

os.makedirs(Fst_path, exist_ok=True)
os.makedirs(Temp_path, exist_ok=True)
os.makedirs(CSV_path, exist_ok=True)
os.makedirs(Means_path, exist_ok=True)

print(f"Created or confirmed: {Fst_path}")
print(f"Created or confirmed: {Temp_path}")
print(f"Created or confirmed: {CSV_path}")
print(f"Created or confirmed: {Means_path}")



# ________________________________________________________________________________________
#                             ###########################
#                                  Jobs in Work-flow  
#                             ###########################
#_________________________________________________________________________________________


#______________________________________STEP 1_____________________________________________

all_FST_logs_vcf = []

# Iterate over individual counts and replication level to define workflow jobs
for replicate in range(no_replicates):
    for no_indv in no_indv_list:
#Draw a defined number of individuals from VCF from two populations.
# + Sort & Zip vcf 
# + Index 
# + Log individual files (metadata): ->  "{Temp_path}/${{pop}}_Tempory_{no_indv}.txt"
        job_name = f"Random_draw_indv_{no_indv}_{replicate}" # Define a unique job name for each individual count and population pair
        job1 = gwf.target_from_template(
            name=job_name + REP_INFO,
            template=Random_hypo_draw(
                VCF=VCF,
                Temp_path = Temp_path + '/' + str(replicate),
                no_indv=no_indv,
                replicate=replicate
            )
        )
#______________________________________STEP 2_____________________________________________

        job2 = gwf.target_from_template(
            name=f"FST_vcftools_{replicate}_{no_indv}" + REP_INFO,#-"-
            template=Calc_FST_vcftools(
                no_miss_hypo_vcf=job1.outputs['no_miss_hypo_vcf'],
                Indvi_txt=job1.outputs['Indvi_txt'],
                Temp_path = Temp_path + '/' + str(replicate),
                no_indv=no_indv, 
                replicate=replicate,
                REP_INFO=REP_INFO,
                Output_path=Fst_path,
                script_path=script_path
            )
        )
        Fst_log=job2.outputs['log']
        all_FST_logs_vcf.append(Fst_log)



# #troubles? - #debug
# #print(','.join(all_FST_logs))


######################################
# #collect the means from .log files.
######################################


#______________________________________STEP 3_____________________________________________

job3 = gwf.target_from_template(
    name='Get_Hypo_FST_mean_vcftools',
    template=Get_hypo_FST_means(
        FSTs=all_FST_logs_vcf,
        Output_path=CSV_path,
        no_indv_list=no_indv_list,
        replicate=replicate,
        Input_path=Fst_path,
        means_path=Means_path
                )
            )

# #troubles? - #debug
# #print(','.join(job3.outputs['Means_file']))

#______________________________________STEP 4_____________________________________________


if clean_workflow == "YES":
    job4 = gwf.target_from_template(
    name='Cleaning',
    template=Clean(
        CSVs=job3.outputs['Means_file'],
        Out_path=Fst_path,
        Temp_path = Temp_path
                    )
                )
else:  # If clean_workflow is not "YES"
    print("Temporary files, VCFs, pr. site Fst, and log files not removed")




