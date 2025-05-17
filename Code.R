# Needed libraries loaded
libraries <- c("dplyr", "skimr", "visdat", "naniar", "mice", "ggplot2", "dbscan", "FactoMineR", "factoextra")
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    library(lib, character.only = TRUE)
  }
}

# Obtain the path of the current working directory
getwd()

#Construct a path relative to the project root
data_path <- "/Users/haoshilong/Downloads/DataSet_No_Details.csv"

# Read the dataset
df <- read.csv(data_path)

# Display the data structure and variable types
str(df)

# Conduct a data overview, including histograms of numeric variables
skim(df) 

# Prepare the dataset
cols_to_remove <- c("h_index_34", "h_index_56", "hormone10_1", "hormone10_2", "an_index_23", "outcome", "factor_eth", "factor_h", "factor_pcos", "factor_prl")
MD_df <- df %>% select(-any_of(cols_to_remove))
factor_df <- df %>% select(record_id, outcome, factor_eth, factor_h, factor_pcos, factor_prl)
str(MD_df)
summary(factor_df)

# Identify missing values
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

# Visualize missing data patterns
vis_miss(MD_df)
gg_miss_var(MD_df)

# Handle missing data by removing some columns
cols_to_remove1 <- c("hormone9", "hormone11", "hormone12", "hormone13", "hormone14")
handle_MD_df <- MD_df %>% select(-any_of(cols_to_remove1))
str(handle_MD_df)

# Perform Little's MCAR test
handle_MD_df_clean <- handle_MD_df %>%
  select(where(~!all(is.na(.)))) %>%
  mutate(across(where(is.character), as.factor))
mcar_result <- mcar_test(handle_MD_df_clean)
print(mcar_result)

# Interpret the MCAR test results
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
# Multiple imputation
# PMM method
imputed_pmm <- mice(handle_MD_df[,!names(handle_MD_df) %in% "New"], method = "pmm")
imputed_pmm_final <- complete(imputed_pmm)

# RF method
imputed_rf <- mice(handle_MD_df[,!names(handle_MD_df) %in% "New"], method = "rf")
imputed_rf_final <- complete(imputed_rf)

# Define a function to compare the density of original data and imputed data
compare_density <- function(var, original_data, pmm_data, rf_data) {
  original <- original_data %>%
    select(all_of(var)) %>%
    mutate(source = "Original")
  pmm <- pmm_data %>%
    select(all_of(var)) %>%
    mutate(source = "PMM")
  rf <- rf_data %>%
    select(all_of(var)) %>%
    mutate(source = "RF")
  combined <- bind_rows(original, pmm, rf)
  
  if (var == "hormone10_generated") {
    x_lim <- c(0, 5)
  } else {
    x_min <- min(combined[[var]], na.rm = TRUE)
    x_max <- max(combined[[var]], na.rm = TRUE)
    x_range <- x_max - x_min
    buffer <- x_range * 0.05
    x_lim <- c(x_min - buffer, x_max - buffer)
  }
  
  ggplot(combined, aes_string(x = var, fill = "source", color = "source")) +
    geom_density(alpha = 0.4, size = 1, na.rm = TRUE) +
    labs(title = paste("Density Comparison of:", var),
         x = var,
         y = "Density") +
    scale_x_continuous(limits = x_lim) +
    scale_fill_manual(values = c(
      "Original" = "black",
      "PMM" = "#E69F00",
      "RF" = "#56B4E9"
    )) +
    scale_color_manual(values = c(
      "Original" = "black",
      "PMM" = "#E69F00",
      "RF" = "#56B4E9"
    )) +
    theme_minimal() +
    theme(legend.position = "top")
}

# Extract the names of variables that need imputation automatically (numeric variables with missing values)
vars_to_plot <- handle_MD_df %>%
  select(where(is.numeric)) %>%
  select(where(~ any(is.na(.)))) %>%
  names()

# Plot density comparison graphs for each variable and save them
for (v in vars_to_plot) {
  plot <- compare_density(v, handle_MD_df, imputed_pmm_final, imputed_rf_final)
  print(plot)
  ggsave(filename = paste0("/Users/ajisafeadediwura/Desktop/Graphs/imputation_density_", v, ".png"), 
         plot = plot, width = 6, height = 4)
}

# Outlier detection
# Select specific columns and reshape the data for outlier detection
library(tidyr)
library(dplyr)
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
ggsave(filename = "/Users/ajisafeadediwura/Desktop/Graphs/outlier_detection_boxplot.png", 
       plot = boxplot1, width = 6, height = 4)

# Generate boxplots for all numeric variables to detect outliers
boxplot2 <- imputed_rf_final %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Boxplots for Outlier Detection")
print(boxplot2)
ggsave(filename = "/Users/ajisafeadediwura/Desktop/Graphs/all_numeric_outlier_boxplot.png", 
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
ggsave(filename = "/Users/ajisafeadediwura/Desktop/Graphs/lof_histogram.png", 
       plot = lof_hist, width = 6, height = 4)

# Scatter plot (show the distribution of LOF scores using lipids1 and lipids2 as an example)
lof_scatter <- ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = lof_score), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "LOF-based Outlier Detection",
       x = "lipids1",
       y = "lipids2",
       color = "LOF Score") +
  theme_minimal()
print(lof_scatter)
ggsave(filename = "/Users/ajisafeadediwura/Desktop/Graphs/lof_scatterplot.png", 
       plot = lof_scatter, width = 6, height = 4)

# Identify outliers (e.g., the top 5% of LOF scores are considered outliers)
threshold <- quantile(lof_scores, 0.95)
lof_df <- lof_df %>%
  mutate(is_outlier = lof_score > threshold)

# Scatter plot highlighting outliers
outlier_scatter <- ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = is_outlier), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(title = "LOF Outliers (Top 5%)",
       color = "Outlier") +
  theme_minimal()
print(outlier_scatter)
ggsave(filename = "/Users/ajisafeadediwura/Desktop/Graphs/lof_outliers_scatterplot.png", 
       plot = outlier_scatter, width = 6, height = 4)
