#negative binomial model to predict outbreak probability

#load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(lmtest)
library(sandwich)
library(zoo)
library(pROC)
library(lme4)
library(DHARMa)

#set working directory if not in project DynamicsofOmicron
setwd('C://Users/srnarum/OneDrive - University of Idaho/CoatsLab/DynamicsofOmicron/Data')
#import data
wastewater <- read_tsv('select_wastewater_data.tsv')
vax <- read_tsv('standardized_vaccines.tsv')
cases <- read_tsv('standardized_cases.tsv')

combined <- full_join(wastewater,cases, by = c("location", "sample_collect_date")) 

#delete unneeded columns
combined <- combined %>%
  select(-daily_cases, -TotalCases)

#create sum of cases that happened in 3 month increments and different vaccine doses in increments
combined <- combined %>%
  group_by(location) %>%  # Group by location
  arrange(sample_collect_date) %>%  # Sort by date
  mutate(
    infected_b = sapply(sample_collect_date, function(date) {
      sum(StandardizedCases[sample_collect_date > (date - days(365)) & sample_collect_date <= (date - days(90))], na.rm = TRUE)
    }),
    infected_c = sapply(sample_collect_date, function(date) {
      sum(StandardizedCases[sample_collect_date > (date - days(455)) & sample_collect_date <= (date - days(365))], na.rm = TRUE)
    }),
    infected_d = sapply(sample_collect_date, function(date) {
      sum(StandardizedCases[sample_collect_date <= (date - days(455))], na.rm = TRUE)
    })
  )

#add vaccines to combined wastewater and cases df
combined.3 <- full_join(combined, vax, by = c("location", "sample_collect_date"), relationship = "many-to-many")
#remove unneeded columns
combined.3 <- combined.3 %>% select(-labels, -County, -Number_of_Vaccines)

#sum vaccines per dose in selected time frames
combined.3 <- combined.3 %>%
  group_by(location) %>%  # Group by location
  arrange(sample_collect_date) %>%  # Sort by date
  mutate(
    two_vax_b = sapply(sample_collect_date, function(date) {
      sum(StandardizedVaccines[sample_collect_date >= (date - days(90)) & 
                                 sample_collect_date < (date - days(14)) & 
                                 DoseNumber == 2], na.rm = TRUE)
    }),
    two_vax_c = sapply(sample_collect_date, function(date) {
      sum(StandardizedVaccines[sample_collect_date >= (date - days(180)) & 
                                 sample_collect_date < (date - days(90)) & 
                                 DoseNumber == 2], na.rm = TRUE)
    }),
    three_vax_b = sapply(sample_collect_date, function(date) {
      sum(StandardizedVaccines[sample_collect_date >= (date - days(60)) & 
                                 sample_collect_date < (date - days(14)) & 
                                 DoseNumber == 3], na.rm = TRUE)
    }),
    three_vax_c = sapply(sample_collect_date, function(date) {
      sum(StandardizedVaccines[sample_collect_date >= (date - days(180)) & 
                                 sample_collect_date < (date - days(60)) & 
                                 DoseNumber == 3], na.rm = TRUE)
    })
  )

#remove rows that are not unique and DoseNumber column
combined.3 <- combined.3 %>% 
  select(-DoseNumber) %>%
  distinct(location, sample_collect_date, .keep_all = TRUE)

#create new columns with protection data from cases and vaccines
#not split by time range
#cases just counting 3-9 months prior
#combined.3$infected_protect <- (combined.3$infected * 0.46)
#vaccines counting 1-6 months ago
#combined.3$boosted_protect <- combined.3$three_vax * 0.672
#combined.3$initial_protect <- combined.3$two_vax * 0.655
#split by time range
#case time range b = 3-12 months post infection
combined.3$infected_b_protect <- (combined.3$infected_b * 0.652) 
#case time range c = 12-15 months post infection
combined.3$infected_c_protect <-(combined.3$infected_c * 0.247) 
#initial time range b = 14 -90 days  post shot
combined.3$two_vax_b_protect <- (combined.3$two_vax_b * 0.44) 
#initial time range c = 90 -180 days post shot
combined.3$two_vax_c_protect <- (combined.3$two_vax_c * 0.235)
#boosted time range b = 14 -60 days  post shot
combined.3$three_vax_b_protect <- (combined.3$three_vax_b * 0.716) 
#boosted time range c = 60 - days post shot
combined.3$three_vax_c_protect <- (combined.3$three_vax_c * 0.474) 


#sum protectiveness amount for different types of events
combined.3$infected_protect <- combined.3$infected_b_protect + combined.3$infected_c_protect
combined.3$boosted_protect <- combined.3$three_vax_b_protect + combined.3$three_vax_c_protect 
combined.3$initial_protect <- combined.3$two_vax_b_protect + combined.3$two_vax_c_protect

#possible option for vax protection calculation
#combined.3$vax_protection <- (combined.3$two_vax_b * 0.66) + (combined.3$two_vax_c * 0.66) + (combined.3$three_vax_b * 0.67) + (combined.3$three_vax_c * 0.67)  
#case protection calculation
#combined.3$case_protection <- (combined.3$infected_b_protect) + (combined.3$infected_c_protect) + (combined.3$infected_d_protect)

#fill in 0 concentration for before omicron came about
combined.3$norm_rolling_mean_omicron[is.na(combined.3$norm_rolling_mean_omicron) & combined.3$sample_collect_date < as.Date("2021-11-01")] <- 0

#remove rows with NA in omicron columns
analysisData <- combined.3 %>%
  subset(!is.na(norm_rolling_mean_omicron)) 

#define outbreak as starting at 3 data points increasing and stopping at 3 data points decreasing

# Sort data by location and date
analysisData <- analysisData %>%
  arrange(location, sample_collect_date)

# Function to determine outbreak
detect_outbreak <- function(concentration) {
  n <- length(concentration)
  outbreak <- rep(0, n)  # Initialize outbreak column with 0s
  
  # Ensure there are no NA values
  concentration <- zoo::na.locf(concentration, na.rm = FALSE)
  
  in_outbreak <- FALSE  # Track if we're in an outbreak
  
  for (i in 3:n) {
    # Detect start of outbreak (2 consecutive increasing values)
    if (!is.na(concentration[i]) && !is.na(concentration[i-1]) && !is.na(concentration[i-2])) {
      if (concentration[i] > concentration[i-1] && concentration[i-1] > concentration[i-2]) {
        in_outbreak <- TRUE  # Outbreak starts
        outbreak[(i-2):i] <- 1
      }
    }
    
    # Mark all samples within outbreak
    if (in_outbreak) {
      outbreak[i] <- 1
    }
    
    # Detect end of outbreak (3 consecutive decreasing values)
    if (i <= (n - 2) && in_outbreak) {
      if (!is.na(concentration[i+1]) && !is.na(concentration[i+2]) && !is.na(concentration[i+3])) {
        if (concentration[i+1] < concentration[i] &&
            concentration[i+2] < concentration[i+1] &&
            concentration[i+3] < concentration[i+2]) {
          in_outbreak <- FALSE  # End outbreak
        }
      }
    }
  }
  
  return(outbreak)
}

# Apply function to each location
analysisData <- analysisData %>%
  group_by(location) %>%
  dplyr::mutate(outbreak = detect_outbreak(norm_rolling_mean_omicron)) %>%
  ungroup()

#Find max omicron concentration values for each location
analysisData_max <- analysisData %>%
  group_by(location) %>%
  dplyr::summarise(max_y = max(norm_rolling_mean_omicron, na.rm = TRUE) + 5, .groups = "drop")
analysisData <- analysisData %>% left_join(analysisData_max, by = "location")

#create ggplot to visualize outbreak samples
ggplot(analysisData, aes(x = sample_collect_date, y = norm_rolling_mean_omicron, color = as.factor(outbreak))) +
  #geom_line() +  # Line plot of rolling mean concentration
  geom_point(size = 2) +  # Points for each sample
  scale_color_manual(values = c("0" = "black", "1" = "red"), name = "Outbreak") +  # Color outbreak points in red
  facet_wrap(~location, scales = "free_y") +  # Separate graphs for each location
  geom_segment(data = analysisData %>% dplyr::filter(outbreak == 1), 
               aes(x = sample_collect_date, xend = sample_collect_date, 
                   y = max_y + 5, yend = max_y + 5), 
               color = "red", size = 2) +  # Horizontal outbreak marker
  labs(title = "Omicron Concentration and Outbreaks", 
       x = "Sample Collection Date", 
       y = "Normalized Rolling Mean Omicron Concentration") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_date(limits = as.Date(c("2021-12-01","2022-02-15")))

#save table
#write.table(analysisData, file = "omicron_outbreak_analysis_data.tsv", col.names = FALSE, row.names = FALSE, sep = "\t")

#remove dates before 07-01-2021 because although protection might be low there is 0% chance of omicron outbreak
analysisData_filtered <- analysisData %>%
  subset(sample_collect_date >= "2021-07-01")
#create  binomial model

#do not have different time ranges for different protectiveness events
#AIC 186
summary(outbreak_with_cases <- glmer(outbreak ~ infected_protect + 
                                boosted_protect + initial_protect + (1|location), 
                              data = analysisData_filtered, 
                              family = binomial))
#AIC 185
summary(outbreak_no_cases <- glmer(outbreak ~ boosted_protect + initial_protect + (1|location), 
                                   data = analysisData_filtered, 
                                   family = binomial))


summary(null_mod <- glmer(outbreak ~ 1 + (1 | location), 
                  data = analysisData_filtered, family = binomial))

#compare models
lrtest(outbreak_with_cases,outbreak_no_cases)

#calculate McFadden's R2
# Extract the log-likelihoods for both models
logLik_full <- logLik(outbreak_no_cases)  # Full model (with predictors)
logLik_null <- logLik(null_mod)    # Null model (intercept-only)

# Calculate McFadden's R²
r2_mcfadden <- 1 - (logLik_full / logLik_null)
r2_mcfadden
#0.67


# Get predicted probabilities
analysisData_filtered$predicted_prob <- predict(outbreak_no_cases, type = "response")
# Convert probabilities to binary predictions
analysisData_filtered$predicted_outbreak <- ifelse(analysisData_filtered$predicted_prob > 0.5, 1, 0)

#auc and roc
roc_curve <- roc(analysisData_filtered$outbreak, analysisData_filtered$predicted_prob)
auc(roc_curve)
#0.984

# Generate the confusion matrix
conf_matrix <- table(Predicted = analysisData_filtered$predicted_outbreak, Actual = analysisData_filtered$outbreak)
print(conf_matrix)
# Compute accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Accuracy:", round(accuracy, 4)))
#0.96
#with cases: 0.96

#compute precision (positive predictive value), recall (sensitivity/true positive rate), F1 score (harmonic mean of precision and recall)
precision <- conf_matrix[2,2] / sum(conf_matrix[2,])
#0.734
#with cases: 0.74
recall <- conf_matrix[2,2] / sum(conf_matrix[,2])
#0.859
#with cases: 0.87
f1_score <- 2 * (precision * recall) / (precision + recall)
#0.79
#with cases:0.80

#residual plots
sim_res <- simulateResiduals(outbreak_no_cases)
plot(sim_res)

#visualize results of the model
#create column  for prediction outcome
analysisData_filtered$prediction_outcome <- with(analysisData_filtered, ifelse(
  outbreak == 1 & predicted_outbreak == 1, "True Positive",
  ifelse(outbreak == 0 & predicted_outbreak == 0, "True Negative",
         ifelse(outbreak == 0 & predicted_outbreak == 1, "False Positive", "False Negative")))
)

#plot the concentration and the true/false negative/positives
plot.prediction <- ggplot(analysisData_filtered, aes(x = sample_collect_date, y = norm_rolling_mean_omicron)) +
  geom_point(aes(color = prediction_outcome), alpha = 1, size = 2.5) +
  geom_line() +
  scale_x_date(limits= as.Date(c("2021-12-01", "2022-02-15"))) +
  facet_wrap(~ location, scales = "free_y") +
  scale_color_manual(values = c(
    "True Positive" = "#1b9e77",
    "True Negative" = "#d95f02",
    "False Positive" = "#7570b3",
    "False Negative" = "#e7298a"
  )) +
  labs(
    x = "Collection Date",
    y = "Normalized Omicron Rolling Mean (copies/day)",
    color = "Prediction Outcome"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
plot.prediction
#ggsave(file=paste("../Figures/Figure4.jpg"), plot=plot.prediction, width=190, height=150, units = "mm", dpi = 300)
