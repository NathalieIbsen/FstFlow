### 
## Visuliaze the estimates


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(tidyr)


##############################################
#1) Exploring the Fst space across populations
##############################################

setwd("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/Cyrtophora_Tammy/")

################# Use "read_and_process_data" to import #################
#
#
complied_FST_indiv_4 <- read_and_process_data("DATA_exploring_fst/Cross_populations/complied_FST_indiv_4.csv", 4)
complied_FST_indiv_6 <- read_and_process_data("DATA_exploring_fst/Cross_populations/complied_FST_indiv_6.csv", 6)
complied_FST_indiv_8 <- read_and_process_data("DATA_exploring_fst/Cross_populations/complied_FST_indiv_8.csv", 8)
complied_FST_indiv_10 <- read_and_process_data("DATA_exploring_fst/Cross_populations/complied_FST_indiv_10.csv", 10)



################# Combine the datasets #################
merged_data <- bind_rows(
  complied_FST_indiv_4,
  complied_FST_indiv_6,
  complied_FST_indiv_8,
  complied_FST_indiv_10
)


################# DATA adjustments #################

# Add a population column to merged_data
merged_data$population <- paste("Cross populations")


# Calculate 
merged_data$no_snps <- cbind(merged_data$no_sites-merged_data$no_monomorph)

################# Plot the create_violin_plot_hypo function #################
# Function to create a violin plot and display mean values as text with several decimals
#You can plot any Fst collum by passing the "collumname" in the function call 

create_violin_plot_hypo(merged_data, "Weir_and_Cockerham_Fst")

create_boxplot_hypo(merged_data, "Weir_and_Cockerham_Fst")


#See how many replicates have above 1000 SNPs
count_reps_above_threshold(merged_data, 2)
count_reps_above_threshold(merged_data, 3)
count_reps_above_threshold(merged_data, 4)
count_reps_above_threshold(merged_data, 5)


#plot a histogram of the distribution in no of snps of each replicate that is going into each mean:

plot1 <- plot_histogram_with_bins(merged_data, 2)
plot2 <- plot_histogram_with_bins(merged_data, 3)
plot3 <- plot_histogram_with_bins(merged_data, 4)
plot4 <- plot_histogram_with_bins(merged_data, 5)


# Adjust relative heights of plots
cowplot::plot_grid(
  plot1, 
  plot2, 
  plot3, 
  plot4, 
  ncol = 5, 
  align = 'v', 
  rel_heights = c(1, 1, 1, 1, 1.3)  # Adjust the relative heights
)


#################Try to remove the run with less that threashold snps, and plot again#################

clean_merged_data <- merged_data %>%
  filter(no_snps > 1000)

create_violin_plot_hypo(clean_merged_data, "Weir_and_Cockerham_Fst")
create_violin_plot_hypo(merged_data, "Weir_and_Cockerham_Fst")

create_boxplot_hypo(clean_merged_data, "Weir_and_Cockerham_Fst")



########################################################
#2) Summarizing and plotting population pairwise Fst 
########################################################


#read in the data

concatenated_runs_all_populations_all_replications <- read.csv("DATA_exploring_fst/All_comp/Min_samples_5_all_comp/concatenated_runs_all_populations_all_replications.csv")

#Calculate a no of SNPs column 

# concatenated_runs_all_populations_all_replications$no_sites <- as.numeric(concatenated_runs_all_populations_all_replications$no_sites)
# 
# concatenated_runs_all_populations_all_replications$no_monomorph <- as.numeric(concatenated_runs_all_populations_all_replications$no_monomorph)


concatenated_runs_all_populations_all_replications$no_snps <- cbind(concatenated_runs_all_populations_all_replications$no_sites-concatenated_runs_all_populations_all_replications$no_monomorph)


# View the result
head(concatenated_runs_all_populations_all_replications)

################# DATA adjustments #################
concatenated_runs_all_populations_all_replications <- concatenated_runs_all_populations_all_replications %>%
  dplyr::rename(no_of_indiv = no_of_indv)


concatenated_runs_all_populations_all_replications <- concatenated_runs_all_populations_all_replications %>%
  mutate(population = ifelse(population == "SPA_CEM", "CEM_SPA", population))

#### Calculate pairwise population means at different sample_sizes 
means_data <- calculate_means(concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")

unique(means_data$population)

table(means_data$population,means_data$no_of_indiv)



##Plot pairwise popualtion means at different sample_sizes
compare_mean_fst(means_data, 2, 3, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(means_data, 3, 4, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(means_data, 4, 5, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(means_data, 2, 5, "mean_Weir_and_Cockerham_Fst")



###plot the distribution of around each mean fst for each pairwise population at all samples sizes
## IF a estimation is done on all samples, no violin is plottet only the "mean". 
#(e.g pop1 has 5 samples, pop2 has 5 samples, if samples size 5 is being compared no violin will appear)


# Usage of the function
create_violin_plot_all(concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")


###But how many SNPs are going in to each estimation?

count_reps_above_threshold(concatenated_runs_all_populations_all_replications, 4)

#Out of how many estimations?  
nrow(concatenated_runs_all_populations_all_replications)

###Try to impose a threshold. 
clean_concatenated_runs_all_populations_all_replications <- concatenated_runs_all_populations_all_replications %>%
  filter(no_snps > 1000)



#How many are left when removing any replicate that didn't have more than 1000 snps 
nrow(clean_concatenated_runs_all_populations_all_replications)

#### Calculate pairwise population means at different sample_sizes 
clean_means_data <- calculate_means(clean_concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")


table(clean_means_data$population,clean_means_data$no_of_indiv)

# View the result
print(clean_means_data)



compare_mean_fst(clean_means_data, 2, 3, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(clean_means_data, 3, 4, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(clean_means_data, 4, 5, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(clean_means_data, 2, 5, "mean_Weir_and_Cockerham_Fst")
compare_mean_fst(clean_means_data, 3, 5, "mean_Weir_and_Cockerham_Fst")


##############################
# For writing out data files.  
# ##############################

# # Splitting the `population` column
# clean_means_data <- clean_means_data %>%
#   separate(population, into = c("pop1", "pop2"), sep = "_")
# 
# 
# # Filter out rows where regions is 'x', 'y', or 'z'
# subset_clean_means_data <- clean_means_data %>%
#   filter(!pop2 %in% c("LAG", "PAL")) %>%
#   filter(!pop1 %in% c("LAG", "PAL"))
# 
# subset_clean_means_data <- subset_clean_means_data %>%
#   filter(!no_of_indiv %in% 3)
# 
# table(subset_clean_means_data$pop1)
# 
# print(subset_clean_means_data)
# 
# #Write to a file: 
# 
# write.csv(subset_clean_means_data, "Sam_mean_Fst_samlesize_5_minSNP200.csv", row.names = FALSE)
# 
# 
# Trial <- clean_concatenated_runs_all_populations_all_replications %>%
#   filter(population == "CCU_SPA")
# 
# 
# # Splitting the `population` column
# means_data <- means_data %>%
#   separate(population, into = c("pop1", "pop2"), sep = "_")
# 
# 
# # Filter out rows where regions is 'x', 'y', or 'z'
# subset_means_data <- means_data %>%
#   filter(!pop2 %in% c("LAG", "PAL")) %>%
#   filter(!pop1 %in% c("LAG", "PAL"))
# 
# subset_means_data <- subset_means_data %>%
#   filter(!no_of_indiv %in% 2)
# 
# 
# unique(subset_clean_means_data$pop2)


#write.csv(subset_means_data, "Sam_mean_Fst_samlesize_5_minSNP200.csv", row.names = FALSE)


# How are these distributed

Estimate_table_data <- as.data.frame(table(concatenated_runs_all_populations_all_replications$population))

Estimate_table_clean_data <- as.data.frame(table(clean_concatenated_runs_all_populations_all_replications$population))

# Rename the 'freq' column in both data frames
colnames(Estimate_table_data) <- c("pop", "Raw Count")
colnames(Estimate_table_clean_data) <- c("pop", "Clean Count")

# Merge the two data frames on 'VAR1'
combined_data_estimates <- merge(Estimate_table_data, Estimate_table_clean_data, by = "pop", all = TRUE)


combined_data_melted <- melt(combined_data_estimates, id.vars = "pop", 
                             variable.name = "Filtering", 
                             value.name = "Frequency")

ggplot(combined_data_melted, aes(x = pop, y = Frequency, fill = Filtering)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim = c(0, 1000)) +  # Set the visible range for the y-axis
  labs(title = "Unique comparisons before and after removing runs below threshold",
       x = "Population Pair",
       y = "No. of Unique Comparisons") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("skyblue", "yellow4"))

#plot the data again - but without the runs with less than 1000 SNPs (the clean_data)


create_violin_plot_all(clean_concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")

create_boxplot_all(concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")

create_boxplot_all(clean_concatenated_runs_all_populations_all_replications, "Weir_and_Cockerham_Fst")

########################################################
##plot the exploration and estimation in the same plot
## COMPAR to the global mean. 

##OBS the data must be run through the two # DATA adjustments # sections above 


# Combine the dataframes by columns
###only the collums found in the dataframe
final_merged_data <- merged_data[, colnames(concatenated_runs_all_populations_all_replications)]

#merge data
final_combined_data<- rbind(concatenated_runs_all_populations_all_replications, final_merged_data)

#set the plotting order 
final_combined_data$population <- factor(final_combined_data$population,
                                         levels = c("Cross populations", 
                                                    unique(final_combined_data$population[
                                                      final_combined_data$population != "Cross populations"])))


# Create the violin plot
create_boxplot_all(final_combined_data, "Weir_and_Cockerham_Fst")
create_boxplot_all(final_combined_data, "Mean_Negatives_Set_to_Zero")
create_boxplot_all(final_combined_data, "Mean_Threshold_Exclusion")



# Combine the dataframes by columns
###only the collums found in the dataframe
clean_final_merged_data <- clean_merged_data[, colnames(clean_concatenated_runs_all_populations_all_replications)]

#merge data
clean_final_combined_data<- rbind(clean_concatenated_runs_all_populations_all_replications, final_merged_data)

#set the plotting order 
clean_final_combined_data$population <- factor(clean_final_combined_data$population,
                                         levels = c("Cross populations", 
                                                    unique(clean_final_combined_data$population[
                                                      clean_final_combined_data$population != "Cross populations"])))




head(clean_final_combined_data)

# Create the boxplot
create_boxplot_all(clean_final_combined_data, "Weir_and_Cockerham_Fst")
create_boxplot_all(clean_final_combined_data, "Mean_Negatives_Set_to_Zero")
create_boxplot_all(clean_final_combined_data, "Mean_Threshold_Exclusion")



## PLOT Fst results for at certain sample size 

#filter - if needed ? 
min3_means_all_populations <- clean_final_combined_data %>%
  filter(no_of_indiv == 3)

head(min3_means_all_populations)


table(min3_means_all_populations$population)


summary(min3_means_all_populations$Weir_and_Cockerham_Fst)

# Get unique levels from the population column
population_levels <- unique(min3_means_all_populations$population)

# Add "Cross populations" as the first level
population_levels <- c("Cross populations", population_levels[!population_levels %in% "Cross populations"])

# Now set the factor levels
min3_means_all_populations$population <- factor(min3_means_all_populations$population,
                                                levels = population_levels)


#########
# PLOT
#########
create_boxplot_all(min3_means_all_populations, "Weir_and_Cockerham_Fst")






#### create a plot with distance on the x 

##First input distance:

#get distance in::: 
# Step 1: Read in the CSV data (adjust paths accordingly)
pairwise_distances <- read.csv("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/Cyrtophora_Tammy/Meta_data/DataPopDistance.csv", sep = ';')

# Step 2: Inspect the data (optional)
head(pairwise_distances)



# Extract unique population pairs from min2_means_all_populations
unique_pairs <- unique(clean_concatenated_runs_all_populations_all_replications$population)
# Create a data frame with unique_pairs as a single column
unique_pairs_df <- data.frame(
  unique_pairs = unique_pairs
)
# Split the unique_pairs into Population1 and Population2 for lookup purposes
unique_pairs_df$Population1 <- sapply(strsplit(unique_pairs, "_"), `[`, 1)
unique_pairs_df$Population2 <- sapply(strsplit(unique_pairs, "_"), `[`, 2)
# Add a column for distances by looking up values in the pairwise_distances matrix
unique_pairs_df$Distance <- mapply(function(p1, p2) {
  # Extract distance for each pair from the matrix
  pairwise_distances[pairwise_distances$X == p1, p2]
}, unique_pairs_df$Population1, unique_pairs_df$Population2)
# Keep only the unique_pairs and Distance columns
final_df <- unique_pairs_df[, c("unique_pairs", "Distance")]

###
#Check on the filtered data#
############################

# View the final DataFrame
head(final_df)
dim(final_df)

dim(clean_concatenated_runs_all_populations_all_replications)
head(clean_concatenated_runs_all_populations_all_replications)


# Merge the data frames using the 'population' column as the key

clean_merged_df <- merge(clean_final_combined_data, final_df, by.x = "population", by.y = "unique_pairs", all.x = TRUE)
clean_merged_df[is.na(clean_merged_df)] <- 0
table(clean_merged_df$population)

merged_df_min3 <- merge(min3_means_all_populations, final_df, by.x = "population", by.y = "unique_pairs", all.x = TRUE)
merged_df_min3[is.na(merged_df_min3)] <- 0




###
#Check on the merged data#
##########################

# View the first few rows of the merged data frame
head(merged_df_min3)
# Check the dimensions of the merged data frame
dim(merged_df_min3)

table(clean_merged_df$population)

################
####NOW PLOT!!!!
################

create_violin_by_factor_distance(clean_merged_df, "Weir_and_Cockerham_Fst")
create_violin_by_factor_distance(merged_df_min3, "Weir_and_Cockerham_Fst")


create_boxplot_by_factor_distance(clean_merged_df, "Weir_and_Cockerham_Fst")
create_boxplot_by_factor_distance(merged_df_min3, "Weir_and_Cockerham_Fst")


create_boxplot_by_numeric_distance(merged_df_min3, "Weir_and_Cockerham_Fst")
create_boxplot_by_numeric_distance(clean_merged_df, "Weir_and_Cockerham_Fst")



#Remove certain populations

# # Splitting the `population` column
# means_data <- means_data %>%
#   separate(population, into = c("pop1", "pop2"), sep = "_")
# 
# 
# # Filter out rows where regions is 'x', 'y', or 'z'
# subset_means_data <- means_data %>%
#   filter(!pop2 %in% c("LAG", "PAL")) %>%
#   filter(!pop1 %in% c("LAG", "PAL"))




















###############
#Subsetting the data 
# aim to color by within region and cross regions. 

DataPopInfo <- read_delim("Meta_data/DataPopInfo.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)

head(final_df)


# Separate unique_pairs into two columns: Pop1 and Pop2
final_df <- final_df %>%
  separate(unique_pairs, into = c("Pop1", "Pop2"), sep = "_")

# Merge to get the side information for Pop1 and Pop2
regional_df <- final_df %>%
  left_join(DataPopInfo %>% select(Pop, Side), by = c("Pop1" = "Pop")) %>%
  rename(Side1 = Side) %>%
  left_join(DataPopInfo %>% select(Pop, Side), by = c("Pop2" = "Pop")) %>%
  rename(Side2 = Side) %>%
  # Combine Side1 and Side2 into the regions column
  mutate(regions = paste(Side1, Side2, sep = "_")) %>%
  # Select the final columns
  mutate(population = paste(Pop1, Pop2, sep = "_")) %>%
  select(population, Distance, regions)

# View the resulting dataframe
head(regional_df)
dim(regional_df)

unique(regional_df$regions)

## as some unique values are the same, eg East_West = West_East I will remap over the names. 

# Transforming the scores
regional_df <- regional_df %>%
  mutate(
    regions = ifelse(regions == "LAG_East", "East_LAG", ifelse(regions == "West_East", "East_West", ifelse(regions == "LAG_West", "West_LAG", ifelse(regions == "PAL_East", "East_LAG", regions)))))

final_merged_df <- merge(concatenated_runs_all_populations_all_replications, regional_df, by.x = "population", by.y = "population", all.x = TRUE)

unique(regional_df$regions)

table(regional_df$regions)

# Merge the data frames using the 'population' column as the key
final_clean_merged_df <- merge(concatenated_runs_all_populations_all_replications, regional_df, by.x = "population", by.y = "population", all.x = TRUE)


# Define a custom color palette
region_colors <- c(
  "East_LAG" = "yellow2",  
  "East_East" = "yellow3", 
  "East_West" = "#2ca02c", 
  "West_West" = "#17fecf",
  "West_LAG" = "#19feef", 
  "West_PAL" = "#19cecf", 
  "East_PAL" = "yellow3", 
  "PAL_LAG" = "darkgray"   
) 

ggplot(final_merged_df, aes(x = Distance, y = Weir_and_Cockerham_Fst, fill = factor(regions))) +
  geom_violin(aes(group = Distance), trim = FALSE, alpha = 0.3, color = "black", 
              adjust = 2, width = 0.8) +  # Adjust bandwidth for smoother violins
  # Plotting smaller dots at the mean with a centered position
  stat_summary(aes(color = factor(regions)), 
               fun = "mean", geom = "point", 
               size = 2) +  # Smaller dots
  stat_summary(aes(color = factor(regions)), 
               fun = function(x) quantile(x, 0.25), geom = "point", 
               shape = 1, size = 2) +
  stat_summary(aes(color = factor(regions)), 
               fun = function(x) quantile(x, 0.75), geom = "point", 
               shape = 1, size = 2) +
  labs(title = "Violin Plot of Weir_and_Cockerham_Fst by Distance",
       x = "Distance",
       y = "Weir_and_Cockerham_Fst",
       fill = "Populations within regions compared",
       color = "Populations within regions compared") +  # Legend for color
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                               margin = margin(t = 4)),  # Smaller x-axis labels
    axis.title.x = element_text(size = 8),  # Smaller x-axis title
    axis.title.y = element_text(size = 8),  # Smaller y-axis title
    legend.text = element_text(size = 6),  # Smaller legend text
    legend.title = element_text(size = 6),  # Smaller legend title
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal",  # Align legend items horizontally
    legend.key.size = unit(0.8, "lines")  # Adjust the size of the legend keys
  ) +
  scale_x_continuous(expand = c(0.01, 0.01)) +  # Numeric x-axis
  ylim(-0.2, 0.6) +
  scale_fill_manual(values = region_colors) +  # Custom fill colors
  scale_color_manual(values = region_colors)  # Custom outline colors



# Create the violin plot with numeric distance
  ggplot(final_merged_df, aes(x = Distance, y = Weir_and_Cockerham_Fst, fill = factor(regions))) +
    geom_violin(aes(group = Distance), trim = FALSE, alpha = 0.3, color = "black", 
                adjust = 2, width = 0.8) +  # Adjust bandwidth for smoother violins
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(regions)), 
                 fun = "mean", geom = "point", 
                 size = 2) +  # Smaller dots
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2) +
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2) +
    labs(title = paste("Violin Plot of Weir_and_Cockerham_Fst by Distance"),
         x = "Distance",
         y = "Weir_and_Cockerham_Fst",
         fill = "Populations within regions compared",
         color = "Populations within regions compared") +  # Legend for color
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  # Smaller x-axis labels
      axis.title.x = element_text(size = 8),  # Smaller x-axis title
      axis.title.y = element_text(size = 8),  # Smaller y-axis title
      legend.text = element_text(size = 6),  # Smaller legend text
      legend.title = element_text(size = 6),  # Smaller legend title
      legend.position = "bottom",  # Move legend to the bottom
      legend.box = "horizontal",  # Align legend items horizontally
      legend.key.size = unit(0.8, "lines")  # Adjust the size of the legend keys
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +  # Numeric x-axis
    ylim(-0.2, 0.6)

  
  

  # Filter out rows where regions is 'x', 'y', or 'z'
subset_final_merged_data <- final_merged_df %>%
    filter(!regions %in% c("East_PAL", "East_LAG","PAL_LAG", "PAL_East", "LAG_East", "West_LAG", "West_PAL"))

subset_region_colors <- c(
  "East_East" = "yellow3", 
  "East_West" = "#2ca02c", 
  "West_West" = "#17fecf")
  
  
  # Create the violin plot with numeric distance
  ggplot(subset_final_merged_data, aes(x = Distance, y = Weir_and_Cockerham_Fst, fill = factor(regions))) +
    geom_violin(aes(group = Distance), trim = FALSE, alpha = 0.3, color = "black", 
                adjust = 2, width = 0.8) +  # Adjust bandwidth for smoother violins
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(regions)), 
                 fun = "mean", geom = "point", 
                 size = 2) +  # Smaller dots
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2) +
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2) +
    labs(title = paste("Violin Plot of Weir_and_Cockerham_Fst by Distance"),
         x = "Distance",
         y = "Weir_and_Cockerham_Fst",
         fill = "Populations within regions compared",
         color = "Populations within regions compared") +  # Legend for color
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  # Smaller x-axis labels
      axis.title.x = element_text(size = 8),  # Smaller x-axis title
      axis.title.y = element_text(size = 8),  # Smaller y-axis title
      legend.text = element_text(size = 6),  # Smaller legend text
      legend.title = element_text(size = 6),  # Smaller legend title
      legend.position = "bottom",  # Move legend to the bottom
      legend.box = "horizontal",  # Align legend items horizontally
      legend.key.size = unit(0.8, "lines")  # Adjust the size of the legend keys
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +  # Numeric x-axis
    ylim(-0.2, 0.6) +
    scale_x_discrete(expand = c(0.01, 0.01)) +  # Reduce space between points
    ylim(-0.15, 0.6) +
    scale_x_continuous(expand = c(0.01, 0.01)) +  # Numeric x-axis
    ylim(-0.2, 0.6) +
    scale_fill_manual(values = subset_region_colors) +  # Custom fill colors
    scale_color_manual(values = subset_region_colors)  # Custom outline colors


  
  
  
  # Create the violin plot with adjusted bandwidth
  ggplot(subset_final_merged_data, aes_string(x = "factor(Distance)", y = "Weir_and_Cockerham_Fst", fill = "factor(regions)")) +
    geom_violin(trim = FALSE, alpha = 0.3, color = "black", 
                position = position_dodge(0.8), width = 1.2, adjust = 2) +  # Adjust bandwidth for smoother violins
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(regions)), 
                 fun = "mean", geom = "point", 
                 size = 2,  # Smaller dots
                 position = position_dodge(0.8)) + # Keep dots centered by using position_dodge
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.05), geom = "point", 
                 shape = 9, size = 2, 
                 position = position_dodge(0.8)) +
    stat_summary(aes(color = factor(regions)), 
                 fun = function(x) quantile(x, 0.95), geom = "point", 
                 shape = 10, size = 2, 
                 position = position_dodge(0.8)) +
    labs(title = paste("Violin Plot ofWeir_and_Cockerham_Fstby Distance"),
         x = "Distance",
         y = "Weir_and_Cockerham_Fst",
         fill = "No. of Individuals",
         color = "No. of Individuals") +  # Legend for color
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  # Smaller x-axis labels
      axis.title.x = element_text(size = 8),  # Smaller x-axis title
      axis.title.y = element_text(size = 8),  # Smaller y-axis title
      legend.text = element_text(size = 6),  # Smaller legend text
      legend.title = element_text(size = 6),  # Smaller legend title
      legend.position = "bottom",  # Move legend to the bottom
      legend.box = "horizontal",  # Align legend items horizontally
      legend.key.size = unit(0.8, "lines")  # Adjust the size of the legend keys
    ) 
  
  