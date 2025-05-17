### Data processing: It is a crucial step in data analysis that involves cleaning, transforming, and organizing raw data into structured format suitable for analysis.
### Data cleaning: The main task is to recognise, identify and handle issues such as errors, missing values, duplicate values, and outliers in the data, ensuring the accuracy, completeness and consistency of the data.

# Introduction:
## 1. Loading necessary libraries
````
libraries <- c("dplyr", "skimr", "visdat", "naniar", "mice", "ggplot2", "dbscan", "FactoMineR", "factoextra")
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    library(lib, character.only = TRUE)
  }
}
````
### The following code defines a set of required R packages and ensures they are installed and loaded. If any package is missing, it is installed and loaded automatically. This setup enables efficient data manipulation, visualization, and analysis.

## 2. Data reading and basic information viewing
```
data_path <- "DataSet_No_Details.csv"
df <- read.csv(data_path)
str(df)
skim(df) 
```
### Specify the dataset file path and read the CSV file into the data frame df. str(df) is used to view the structure of the data frame, understand the data types and general content of each column; skim(df) provides a detailed summary of the data, such as the distribution of numerical variables.

## 3. Dataset preparation
```
cols_to_remove <- c("h_index_34", "h_index_56", "hormone10_1", "hormone10_2", "an_index_23", "outcome", "factor_eth", "factor_h", "factor_pcos", "factor_prl")
MD_df <- df %>% select(-any_of(cols_to_remove))
factor_df <- df %>% select(record_id, outcome, factor_eth, factor_h, factor_pcos, factor_prl)
str(MD_df)
summary(factor_df)
```
### Specify a vector of column names to exclude from the original data frame df, and utilize the dplyr package's select function to create MD_df by removing these columns. Additionally, extract specific columns to form factor_df. Finally, examine the structure of MD_df and generate summary statistics for factor_df.

## 4. Missing value identification and analysis
```
total_na <- sum(is.na(MD_df))               
col_na_counts <- colSums(is.na(MD_df))           
skim(MD_df)
na_stats <- colMeans(is.na(MD_df)) * 100 
na_stats_filtered <- na_stats[na_stats <= 35] 
na_stats_filtered_table <- data.frame(
  Column = names(na_stats_filtered),
  NA_Percent = na_stats_filtered,
  row.names = NULL
)
na_stats_filtered_1 <- na_stats[na_stats > 35] 
na_stats_filtered_1_table <- data.frame(
  Column = names(na_stats_filtered_1),
  NA_Percent = na_stats_filtered_1,
  row.names = NULL
)
vis_miss(MD_df)
gg_miss_var(MD_df)
```
### Compute the total count of missing values in MD_df, as well as the count and percentage of missing values for each column. Then, categorize columns based on their missing value percentages: one group with percentages not exceeding 35% and another with percentages above 35%. Organize these columns into separate data frames. Finally, utilize the vis_miss and gg_miss_var functions to visualize the patterns of missing data.

## 5. Missing data processing and MCAR test
```
cols_to_remove1 <- c("hormone9", "hormone11", "hormone12", "hormone13", "hormone14")
handle_MD_df <- MD_df %>% select(-any_of(cols_to_remove1))
handle_MD_df_clean <- handle_MD_df %>%
  select(where(~!all(is.na(.)))) %>%
  mutate(across(where(is.character), as.factor))
mcar_result <- mcar_test(handle_MD_df_clean)
print(mcar_result)
interpret_mcar <- function(mcar_result) {
  p <- mcar_result$p.value
  if (p > 0.05) {
    message("✅ p-value > 0.05 → Data is likely MCAR. Safe to delete or impute.")
  } else {
    message("❗ p-value <= 0.05 → Data is NOT MCAR. Assume MAR or MNAR.")
    message("➡️ It is recommended to use the multiple imputation (mice) method, such as pmm / rf, etc.")
  }
}
interpret_mcar(mcar_result)
```
### Further refine MD_df by removing columns with excessive missing values to create handle_MD_df. Then, process handle_MD_df to obtain handle_MD_df_clean by eliminating columns with all missing values and converting character variables to factor variables. Perform Little's MCAR test to determine if the data is Missing Completely At Random (MCAR). Based on the test results, assess whether the data meets the MCAR assumption and provide recommendations for subsequent data processing.

## 6. Multiple imputation
```
imputed_pmm <- mice(handle_MD_df[,!names(handle_MD_df) %in% "New"], method = "pmm")
imputed_pmm_final <- complete(imputed_pmm)
imputed_rf <- mice(handle_MD_df[,!names(handle_MD_df) %in% "New"], method = "rf")
imputed_rf_final <- complete(imputed_rf)
```
### Utilize the mice package to perform multiple imputations on the processed data using two methods: Predictive Mean Matching (PMM) and Random Forest (RF). Generate complete datasets with imputed values using both approaches.
## 7. Comparison of imputation effects
```
compare_density <- function(var, original_data, pmm_data, rf_data) {
  # Function content omitted
}
vars_to_plot <- handle_MD_df %>%
  select(where(is.numeric)) %>%
  select(where(~ any(is.na(.)))) %>%
  names()
for (v in vars_to_plot) {
  plot <- compare_density(v, handle_MD_df, imputed_pmm_final, imputed_rf_final)
  print(plot)
  ggsave(filename = paste0("Users/haoshilong/Desktop/Graphs/imputation_density_", v, ".png"), 
         plot = plot, width = 6, height = 4)
}
```
### Create a function compare_density to visualize and compare the density distributions of the original data with those obtained from two imputation methods. Identify numerical variables with missing values, generate density plots for each variable, and save the plots to a specified directory.

## 8. Outlier detection
```# Outlier detection
# Select specific columns and reshape the data for outlier detection
outliers_data <- imputed_rf_final %>%
  select(lipids1, lipids2, lipids3, lipids4, lipids5) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")
# Generate boxplots to detect outliers
boxplot1 <- ggplot(outliers_data, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Outlier Detection",
       x = "variables",
       y = "value") +
  theme_minimal()
print(boxplot1)
ggsave(filename = "Users/haoshilong/Desktop/Graphs/outlier_detection_boxplot.png", 
       plot = boxplot1, width = 6, height = 4)
# Generate boxplots for all numerical variables to detect outliers
boxplot2 <- imputed_rf_final %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Boxplots for Outlier Detection")
print(boxplot2)
ggsave(filename = "Users/haoshilong/Desktop/Graphs/all_numeric_outlier_boxplot.png", 
       plot = boxplot2, width = 8, height = 6)
# Calculate LOF factors
lof_data <- imputed_rf_final %>%
  select(where(is.numeric)) %>%
  na.omit()
lof_scores <- lof(lof_data, k = 20)
lof_df <- lof_data %>%
  mutate(lof_score = lof_scores)
# Visualize LOF factors
# Histogram
lof_hist <- ggplot(lof_df, aes(x = lof_score)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "#FF7F00",
    color = "black",
    alpha = 0.7
  ) +
  labs(
    title = "Histogram of  LOF Scores",
    x = "LOF Score",
    y = "Frequency"
  ) +
  theme_minimal()
print(lof_hist)
ggsave(filename = "Users/haoshilong/Desktop/Graphs/lof_histogram.png", 
       plot = lof_hist, width = 6, height = 4)
# Scatter plot (show the distribution of LOF scores using lipids1 and lipids2 as an example)
lof_scatter <- ggplot(lof
```
### Analyze the imputed data for outliers using box plots and the Local Outlier Factor (LOF) algorithm. Generate various box plots to visualize data distributions, calculate LOF scores to identify potential outliers, and create histograms and scatter plots to visualize these scores. Determine outliers based on LOF scores, highlight them in the plots, and save all relevant graphics to the specified directory.
