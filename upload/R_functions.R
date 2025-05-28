### R function to handle the output of:


########################################################
#1) Exploring the Fst space across populations
########################################################

########################################################
#2) Summarizing and plotting population pairwise Fst 
########################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)

###### 1 FUNCTION ######

# Read in the data 

read_and_process_data <- function(file_path, no_of_indiv) {
  data <- read.csv(file_path)
  
  data <- data %>%
    mutate(no_of_indiv = no_of_indiv / 2)  # Divide the `no_of_indiv` by 2
  data$no_snps <- cbind(as.numeric(data$no_sites-data$no_monomorph))
  
  return(data)
}


###### 2 FUNCTION ######


# Function to create a violin plot and display mean values as text with several decimals
create_violin_plot_hypo <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Ensure no_of_indiv is treated as a factor to get one violin plot for each unique value
  data$no_of_indiv <- as.factor(data$no_of_indiv)  # Correct column name
  
  # Create the violin plot
  plot <- ggplot(data, aes(x = no_of_indiv, y = !!sym(column_to_plot))) +
    geom_violin(trim = FALSE, fill = "lightblue") +  # Violin plot
    stat_summary(fun = mean, geom = "point", 
                 color = "blue", size = 3, 
                 position = position_dodge(0.9), 
                 na.rm = TRUE) +  # Ensures that missing values are ignored
    labs(
         x = "Number of Individuals",
         y = column_to_plot) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Optional theme for better aesthetics
  
  # Add mean value as text with several decimals
  plot + geom_text(
    data = data %>%
      group_by(no_of_indiv) %>%
      summarise(mean_value = mean(!!sym(column_to_plot), na.rm = TRUE)),  # Calculate the mean for each group
    aes(x = no_of_indiv, y = mean_value, label = sprintf("%.4f", mean_value)),  # Format to 4 decimal places
    color = "blue", 
    vjust = -1, size = 3  # Adjust vertical positioning and text size
  ) +
    ylim(-0.2, 0.25)  # Set y-axis limits
}



###### 2(b) FUNCTION ######


create_boxplot_hypo <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Ensure no_of_indiv is treated as a factor to get one boxplot per unique value
  data$no_of_indiv <- as.factor(data$no_of_indiv)
  
  # Create the boxplot
  plot <- ggplot(data, aes(x = no_of_indiv, y = !!sym(column_to_plot))) +
    geom_boxplot(fill = "lightblue", outlier.color = "lightblue", outlier.shape = 16, outlier.size = 1) +  # Boxplot with red outliers
    stat_summary(fun = mean, geom = "point", 
                 color = "blue", size = 3, 
                 position = position_dodge(0.9), 
                 na.rm = TRUE) +  # Mean points
    labs(
      x = "Number of Individuals",
      y = column_to_plot) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional theme for better aesthetics
  
  # Add mean value as text with several decimals
  plot + geom_text(
    data = data %>%
      group_by(no_of_indiv) %>%
      summarise(mean_value = mean(!!sym(column_to_plot), na.rm = TRUE)),  # Calculate mean per group
    aes(x = no_of_indiv, y = mean_value, label = sprintf("%.4f", mean_value)),  # Format to 4 decimal places
    color = "blue", 
    vjust = -2, size = 2  # Adjust text positioning and size
  ) +
    ylim(-0.2, 0.25)  # Set y-axis limits (adjust as needed)
}





###### 3 FUNCTION ######

count_reps_above_threshold <- function(data, indiv_value, threshold = 999) {
  # Ensure data has expected structure
  if (!"no_of_indiv" %in% names(data) || !"no_snps" %in% names(data)) {
    stop("Input data must contain 'no_of_indiv' and 'no_snps' columns.")
  }
  
  # Filter data for the specified no_of_indiv
  filtered_data <- data %>% filter(no_of_indiv == indiv_value)
  
  # Ensure that the condition creates a logical vector
  count <- filtered_data %>% filter(no_snps > threshold) %>% nrow()
  
  return(count)
}


###### 4 FUNCTION ######


#plot a histogram of the distribution in no of snps of each replicate that is going into each mean:
#in bins of 100

plot_histogram_with_bins <- function(data, indiv_value, bins = c(0, 99, 199, 299, 399, 499, 599, 699, 799, 999, Inf)) {
  # Filter data for the specified no_of_indiv
  filtered_data <- data %>% filter(no_of_indiv == indiv_value)
  
  # Categorize no_snps into the specified bins
  filtered_data <- filtered_data %>%
    mutate(snps_bins = cut(no_snps, breaks = bins, right = TRUE, include.lowest = TRUE))
  
  # Create the histogram
  ggplot(filtered_data, aes(x = snps_bins)) +
    geom_bar(color = "black", fill = "lightblue") +
    labs(
      title = paste("No. of SNPs for", indiv_value, "individuals"),
      x = "SNP Count Bins",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=6, angle = 45, hjust = 1)) +
    ylim(0,900)
}


###### 5 FUNCTION ######


#calculate pairwise population means

calculate_means <- function(data, column_to_average) {
  # Check if the specified column exists in the data frame
  if (!column_to_average %in% names(data)) {
    stop(paste("Column", column_to_average, "not found in the data frame."))
  }
  
  # Calculate the mean for each population and no_of_indiv combination
  result <- data %>%
    group_by(population, no_of_indiv) %>%
    summarise(mean_value = mean(.data[[column_to_average]], na.rm = TRUE), .groups = "drop")
  
  # Rename the mean column for clarity
  result <- result %>%
    rename(!!paste0("mean_", column_to_average) := mean_value)
  
  return(result)
}


###### 5 FUNCTION ######


###Plot sample_no comparison, plot the mean of indv_1 against the mean of indv_2

# Define the function
compare_mean_fst <- function(data, indv_1, indv_2, column) {
  # Convert the column name to a symbol for use in tidy evaluation
  column_sym <- sym(column)
  
  # Filter for the first specified no_of_indiv and rename for merging
  data_1 <- data %>%
    filter(no_of_indiv == indv_1) %>%
    select(population, !!column_sym) %>%
    rename(Value_1 = !!column_sym)
  
  # Filter for the second specified no_of_indiv and rename for merging
  data_2 <- data %>%
    filter(no_of_indiv == indv_2) %>%
    select(population, !!column_sym) %>%
    rename(Value_2 = !!column_sym)
  
  # Join the data on the `population` column
  merged_data <- inner_join(data_1, data_2, by = "population")
  
  # Plot the specified column for the two specified no_of_indiv values
  ggplot(merged_data, aes(x = Value_1, y = Value_2, label = population)) +
    geom_point() +
    geom_text(vjust = -1, hjust = 1, size = 3) +  # Adding population labels to points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Equality line
    labs(
      x = paste(column, "(no_of_indiv =", indv_1, ")"),
      y = paste(column, "(no_of_indiv =", indv_2, ")"),
      title = paste("Comparison of", column, "for no_of_indiv =", indv_1, "vs", indv_2, "by Population")
    ) +
    theme_minimal()
}


###### 6 FUNCTION ######


###plot the distribution of around each mean fst for each pairwise population at all samples sizes  
## IF a estimation is done on all samples, no violin is plottet only the "mean". 
#(e.g pop1 has 5 samples, pop2 has 5 samples, if samples size 5 is being compared no violin will appear)


# Function to create a violin plot
create_violin_plot_all <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Calculate the difference between max and min for each population and no_of_indiv group
  data_summary <- data %>%
    group_by(population, no_of_indiv) %>%
    summarize(diff = max(.data[[column_to_plot]], na.rm = TRUE) - min(.data[[column_to_plot]], na.rm = TRUE),
              .groups = 'drop')
  
  # Add a flag to the data to indicate if it should be plotted as an empty violin
  data <- data %>%
    left_join(data_summary, by = c("population", "no_of_indiv")) %>%
    mutate(plot_violin = ifelse(diff < 0.0001, "empty", "filled"))
  
  # Create the ggplot object
  p <- ggplot(data, aes_string(x = "population", y = column_to_plot, fill = "factor(no_of_indiv)")) +
    # Plot filled violins where the condition is not met
    geom_violin(trim = FALSE, alpha = 0.3, color = "black", 
                position = position_dodge(0.8), width = 1.5, adjust = 1, size = 0.1, 
                data = subset(data, plot_violin == "filled")) +  # Plot filled violins
    
    # For empty violins, plot a blank space (empty violin with no width and no fill)
    geom_blank(data = subset(data, plot_violin == "empty"), 
               aes(x = population, y = 0)) +  # Create empty space for violins
    
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2,  # Smaller dots
                 position = position_dodge(0.8)) + # Keep dots centered by using position_dodge
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    labs(title = paste("Violin Plot of", column_to_plot, "by Population"),
         x = "Population",
         y = gsub("_", " ", column_to_plot),
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
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +  # Reduce space between points
    ylim(-0.1, 0.5)
  
  return(p)
}


###### 6(b) FUNCTION ######



create_boxplot_all <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Calculate the difference between max and min for each population and no_of_indiv group
  data_summary <- data %>%
    group_by(population, no_of_indiv) %>%
    summarize(diff = max(.data[[column_to_plot]], na.rm = TRUE) - min(.data[[column_to_plot]], na.rm = TRUE),
              .groups = 'drop')
  
  # Add a flag to the data to indicate if it should be plotted as an empty boxplot
  data <- data %>%
    left_join(data_summary, by = c("population", "no_of_indiv")) %>%
    mutate(plot_box = ifelse(diff < 0.0001, "empty", "filled"))
  
  # Create the ggplot object
  p <- ggplot(data, aes_string(x = "population", y = column_to_plot, fill = "factor(no_of_indiv)")) +
    # Plot filled boxplots where the condition is met
    geom_boxplot(data = subset(data, plot_box == "filled"), 
                 outlier.shape = NA, alpha = 0.3, color = "black", 
                 position = position_dodge(0.8), width = 0.6) +
    
    # For empty boxplots, plot a blank space (boxplot with no width and no fill)
    geom_blank(data = subset(data, plot_box == "empty"), 
               aes(x = population, y = 0)) +  
    
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2, 
                 position = position_dodge(0.8)) +
    
    labs(title = paste("Boxplot of", column_to_plot, "by Population"),
         x = "Populationpair",
         y = gsub("_", " ", column_to_plot),
         fill = "No. of Individuals",
         color = "No. of Individuals") +  
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  
      axis.title.x = element_text(size = 8),  
      axis.title.y = element_text(size = 8),  
      legend.text = element_text(size = 6),  
      legend.title = element_text(size = 6),  
      legend.position = "bottom",  
      legend.box = "horizontal",  
      legend.key.size = unit(0.8, "lines")  
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    ylim(-0.1, 0.5)
  
  return(p)
}






###### 7 FUNCTION ######
###plot the distribution around each mean fst for each pairwise population, at a certain sample size
## IF a estimation is done on all samples, no violin is plottet only the "mean". 
#(e.g pop1 has 5 samples, pop2 has 5 samples, if samples size 5 is being compared no violin will appear)




# Function to create a violin plot with conditional plotting
create_violin_plot_min <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Group the data by population and check if the range of column_to_plot is sufficiently large
  data_grouped <- data %>%
    dplyr::group_by(population) %>%
    dplyr::mutate(range_column = max(.data[[column_to_plot]], na.rm = TRUE) - 
                    min(.data[[column_to_plot]], na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # Create the base ggplot object
  plot <- ggplot(data_grouped, aes_string(x = "population", y = column_to_plot, fill = "factor(no_of_indiv)")) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2,  # Smaller dots
                 position = position_dodge(0.8)) +  # Keep dots centered by using position_dodge
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    labs(title = paste("Violin Plot of", column_to_plot, "by Population"),
         x = "Population",
         y = gsub("_", " ", column_to_plot),
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
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +  # Reduce space between points
    ylim(-0.1, 0.6)
  
  # Add the violin plot only if the range is greater than the threshold (0.0001)
  plot <- plot + geom_violin(
    data = data_grouped %>% filter(range_column > 0.0001),
    trim = FALSE, alpha = 0.3, color = "black", 
    position = position_dodge(0.8), width = 1.2, adjust = 2)
  
  return(plot)
}



###### 7 FUNCTION ######

###plot the distribution around each mean fst for each pairwise population, BUT in the order of distance, at a certain sample size
## IF a estimation is done on all samples, no violin is plottet only the "mean". 
#(e.g pop1 has 5 samples, pop2 has 5 samples, if samples size 5 is being compared no violin will appear)



# Function to create a violin plot based on Distance
create_violin_by_factor_distance <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Create the violin plot with adjusted bandwidth
  ggplot(data, aes_string(x = "factor(Distance)", y = column_to_plot, fill = "factor(no_of_indiv)")) +
    geom_violin(trim = FALSE, alpha = 0.3, color = "black", 
                position = position_dodge(0.8), width = 1.2, adjust = 2) +  # Adjust bandwidth for smoother violins
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2,  # Smaller dots
                 position = position_dodge(0.8)) + # Keep dots centered by using position_dodge
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2, 
                 position = position_dodge(0.8)) +
    labs(title = paste("Violin Plot of", column_to_plot, "by Distance"),
         x = "Distance",
         y = gsub("_", " ", column_to_plot),
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
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +  # Reduce space between points
    ylim(-0.15, 0.7)
}


###### 7(bS) FUNCTION ######

###plot the distribution around each mean fst for each pairwise population, BUT in the order of distance, at a certain sample size
## IF a estimation is done on all samples, no violin is plottet only the "mean". 
#(e.g pop1 has 5 samples, pop2 has 5 samples, if samples size 5 is being compared no violin will appear)


create_boxplot_by_factor_distance <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Create the boxplot
  ggplot(data, aes_string(x = "factor(Distance)", y = column_to_plot, fill = "factor(no_of_indiv)")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, color = "black", 
                 position = position_dodge(0.8), width = 0.6) +
    
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2,  
                 position = position_dodge(0.8)) + 
    
    labs(title = paste("Boxplot of", gsub("_", " ", column_to_plot), "by Distance"),
         x = "Distance(Km)",
         y = gsub("_", " ", column_to_plot),
         fill = "No. of Individuals",
         color = "No. of Individuals") +  
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  
      axis.title.x = element_text(size = 8),  
      axis.title.y = element_text(size = 8),  
      legend.text = element_text(size = 6),  
      legend.title = element_text(size = 6),  
      legend.position = "bottom",  
      legend.box = "horizontal",  
      legend.key.size = unit(0.8, "lines")  
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    ylim(-0.15, 0.5)
}



create_violin_by_numeric_distance <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Create the violin plot with numeric distance
  ggplot(data, aes(x = Distance, y = !!sym(column_to_plot), fill = factor(no_of_indiv))) +
    geom_violin(aes(group = Distance), trim = FALSE, alpha = 0.3, color = "black", 
                adjust = 2, width = 0.8) +  # Adjust bandwidth for smoother violins
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2) +  # Smaller dots
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2) +
    labs(title = paste("Violin Plot of", column_to_plot, "by Distance"),
         x = "Distance",
         y = column_to_plot,
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
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +  # Numeric x-axis
    ylim(-0.1, 0.5)
}


create_boxplot_by_numeric_distance <- function(data, column_to_plot) {
  # Check if the specified column exists in the dataframe
  if (!column_to_plot %in% colnames(data)) {
    stop(paste("Column", column_to_plot, "does not exist in the dataframe."))
  }
  
  # Create the boxplot with numeric distance
  ggplot(data, aes(x = Distance, y = !!sym(column_to_plot), fill = factor(no_of_indiv))) +
    geom_boxplot(aes(group = Distance), outlier.shape = NA, alpha = 0.3, color = "black", 
                 width = 0.6) +
    
    # Plotting smaller dots at the mean with a centered position
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = "mean", geom = "point", 
                 size = 2) +  
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.25), geom = "point", 
                 shape = 1, size = 2) +
    stat_summary(aes(color = factor(no_of_indiv)), 
                 fun = function(x) quantile(x, 0.75), geom = "point", 
                 shape = 1, size = 2) +
    
    labs(title = paste("Boxplot of", gsub("_", " ", column_to_plot), "by Distance"),
         x = "Distance",
         y = gsub("_", " ", column_to_plot),
         fill = "No. of Individuals",
         color = "No. of Individuals") +  
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, 
                                 margin = margin(t = 4)),  
      axis.title.x = element_text(size = 8),  
      axis.title.y = element_text(size = 8),  
      legend.text = element_text(size = 6),  
      legend.title = element_text(size = 6),  
      legend.position = "bottom",  
      legend.box = "horizontal",  
      legend.key.size = unit(0.8, "lines")  
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    ylim(-0.1, 0.6)
}


