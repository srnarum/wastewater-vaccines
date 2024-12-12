#load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(lmtest)
library(psc1)

#set working directory if not in project DynamicsofOmicron
setwd('C://Users/srnarum/OneDrive - University of Idaho/CoatsLab/DynamicsofOmicron/Data')
#import data
data <- read_tsv('combined_data.tsv')
earliest_dates <- read_tsv('omicron_earliest_detection.tsv')
cases <- read_tsv('all_cases_data.tsv')

#change units of concentration of delta, omicron, N1 to be (gc/ml) instead of gc/uL
data$rolling_mean_delta <- data$rolling_mean_delta * 1000
data$rolling_mean_omicron <- data$rolling_mean_omicron * 1000
data$rolling_mean_N <- data$rolling_mean_N * 1000

#create df with cumulative sum of vaccines standardized to location per shot
percent.vax <- data %>% 
  mutate(shot1_percent = (shot_1 * 100) / population) %>%
  mutate(shot2_percent = (shot_2 * 100) / population) %>%
  mutate(shot3_percent = (shot_3 * 100) / population) %>%
  mutate(shot4_percent = (shot_4 * 100) / population) 

#order by location and date
vaccinated <- percent.vax[order(percent.vax$location, percent.vax$sample_collect_date), ]
# Unique locations
locations <- unique(vaccinated$location)
# Initialize the initialShots column
vaccinated$initialShots <- NA

# Calculate sum of shot 3 shots and subtract after 180 days for each row
# Calculate boosted1 for each location separately
for (loc in locations) {
  loc_indices <- which(vaccinated$location == loc)
  loc_data <- vaccinated[loc_indices, ]
  
  for (i in seq_along(loc_indices)) {
    current_date <- loc_data$sample_collect_date[i]
    start_date <- current_date - 180
    relevant_rows <- loc_data$sample_collect_date <= current_date & loc_data$sample_collect_date > start_date
    vaccinated$initialShots[loc_indices[i]] <- sum(loc_data$shot2_percent[relevant_rows])
  }
}

# Initialize the boosted1 column
vaccinated$boosted1 <- NA

# Calculate sum of shot 3 shots and subtract after 180 days for each row
# Calculate boosted1 for each location separately
for (loc in locations) {
  loc_indices <- which(vaccinated$location == loc)
  loc_data <- vaccinated[loc_indices, ]
  
  for (i in seq_along(loc_indices)) {
    current_date <- loc_data$sample_collect_date[i]
    start_date <- current_date - 180
    relevant_rows <- loc_data$sample_collect_date <= current_date & loc_data$sample_collect_date > start_date
    vaccinated$boosted1[loc_indices[i]] <- sum(loc_data$shot3_percent[relevant_rows])
  }
}

#create up_to_date column which includes initial round of shots and vaccinated for plotting
vaccinated$up_to_date <- vaccinated$initialShots + vaccinated$boosted1

# Now plot the data
up_to_date.plot <- ggplot(vaccinated, aes(x = sample_collect_date, y = up_to_date)) +
  geom_point() +
  geom_vline(data = earliest_dates, aes(xintercept = as.numeric(earliest_date)), color = "blue")+
  geom_vline(aes(xintercept = as.Date('2021-10-31')), color = "orange")+
  geom_vline(aes(xintercept = as.Date('2022-02-15')), color = 'orange')+
  labs(x = "Date", y = "Percent of Population", title = "Protection Against Omicron (including initial vaccine and booster)") +
  facet_wrap(~ location)+
  scale_x_date(date_breaks = "1 month",
               limits = as.Date(c("2021-07-01", "2022-07-01"))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10))

up_to_date.plot

#save the plot
#ggsave(filename = '../Figures/vaccine_protection_from_omicron.jpg', plot = up_to_date.plot)

#prepare case data with protectiveness from infection data
#take infections that occured within 13 months before Omicron start and multiply by 36% to get immunity to omicron

#just plot one town
Moscow <- vaccinated %>%
  subset(location == "Moscow")

Moscow.plot <- ggplot(Moscow, aes(x = sample_collect_date, y = up_to_date)) +
  geom_point(color = "purple") +
  geom_vline(aes(xintercept = as.numeric("2021-12-12")), color = "blue")+
  labs(x = "", y = "") +
  scale_x_date(date_breaks = "1 month",
               limits = as.Date(c("2021-07-01", "2022-07-01")),
               date_labels = "%b %Y") +
 theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank())
Moscow.plot
#ggsave(filename = '../Figures/Moscow_protection_from_omicron.jpg', plot = Moscow.plot, width = 3, height = 1.5)

#remove rows with NA (only want dates with wastewater data)
analysisData <- na.omit(vaccinated)

#3 parameters
summary(mod.3 <- glm(rolling_mean_omicron ~ initialShots * boosted1 + offset(log(rolling_mean_N)), family = gaussian(link = "log"), 
            data = analysisData, mustart = pmax(analysisData$rolling_mean_omicron,1e-3)))
#2 parameters
#this shows similar residual and null deviances as with * , but also has signigicant p-values for the parameters
summary(mod.2 <- glm(rolling_mean_omicron ~ initialShots + boosted1 + offset(log(rolling_mean_N)), family = gaussian(link = "log"), 
            data = analysisData, mustart = pmax(analysisData$rolling_mean_omicron,1e-3)))
#0 parameters
summary(mod.0 <- glm(rolling_mean_omicron ~ 1 + offset(log(rolling_mean_N)), family = gaussian(link = "log"), 
            data = analysisData, mustart = pmax(analysisData$rolling_mean_omicron,1e-3)))

#to compare models
lrtest(mod.0, mod.2)
# 2.2e-16

lrtest(mod.2, mod.3)
#0.7716

#determine how much better 2 parameter model is to 0 parameter model
AIC(mod.2,mod.0)
BIC(mod.2,mod.0)
#discuss 3 possible models, compare with significance using likelihood ratio test, 
#give pseudo-R squared and significance of chisq test and table of coefficients and p values for the best model

#plot predictions of model vs observed data
# Calculate the predictions
analysisData$predictions <- predict(mod.2, type="response")

# Define ylims for matching with the base R code's plot limits
plot(analysisData$rolling_mean_omicron, predict(mod.2, type="response"), log = "xy")
ylims <- 10^par("usr")[1:2]

#add a small constant so there are no zero values
analysisData$rolling_mean_omicron <- analysisData$rolling_mean_omicron + 1e-6
analysisData$predictions <- analysisData$predictions + 1e-6
#remove data with small observed values
analysisData_filtered <- analysisData[analysisData$rolling_mean_omicron > 0.000002, ]
# Create the plot using ggplot2
pred.vs.obs <- ggplot(analysisData_filtered, aes(x = rolling_mean_omicron, y = predictions)) +
  geom_point() +  # Scatter plot for the data points
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  labs(
    x = "log(Observed Omicron \n Concentration (gc/ml))",
    y = "log(Prediction of Omicron \n Concentration (gc/ml))"
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Red identity line
  annotate("text", 
           x = 50, 
           y = 30, 
           label = expression(italic(rho) == 0.81 * ", " * hat(y) * " = " * - 0.165 - 0.471 ~ x[1] + 0.108 ~ x[2]), 
           size = 3, 
           hjust = 0, 
           vjust = 1) +  # Annotating the formula
  #coord_cartesian(ylim = ylims) +  # Set y-axis limits based on ylims
  #theme_minimal(base_family = "Arial") + # Setting the font family to Arial
  theme(panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", size = 1),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7))

pred.vs.obs
ggsave("../Figures/Figure4.jpg", plot = pred.vs.obs, height = 70, width = 90, dpi = 300, units = "mm")

#plot the correlation between the percentage of the population protected against omicron due to the initial vaccine 
#and the concentration of omicron in the wastewater
plot.correlation = ggscatter(data = analysisData_filtered, x = "initialShots", y = "rolling_mean_omicron", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "Protection due to Vaccine", ylab = "7-day Rolling Mean of Omicron Concentration (gc/ml)")

plot.correlation

#plot the correlation between the percentage of the population that received the booster shot 
#and the concentration of omicron in the wastewater
plot.correlation = ggscatter(data = analysisData_filtered, x = "boosted1", y = "rolling_mean_omicron", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "Protection due to Booster", ylab = "7-day Rolling Mean of Omicron Concentration (gc/ml)")

plot.correlation

#previous method of plotting predicted vs observed, but in base R
plot(analysisData$rolling_mean_omicron, predict(mod.2, type="response"), log = "xy")
ylims <- 10^par("usr")[1:2]

#jpeg("../Pub/Environmental Science and Technology/Figure4.jpg", res = 300, width = 4, height = 3, units = "in", pointsize = 8)
#adjust margins
par(mar = c(5, 5, 0.5, 0.5), family = "Arial")
plot(analysisData$rolling_mean_omicron, predict(mod.3, type="response"), log="xy", 
     xlab="log(Observed Omicron Concentration (gc/ml))",
     ylab = "log(Prediction of Omicron Concentration(gc/ml))", ylim = ylims)
abline(a=0,b=1, col="red")
text(0.1,5, expression(italic(rho) == 0.81 * ", " * hat(y) * " = " * - 0.165 - 0.471 ~ x[1] + 0.108 ~ x[2]))
#dev.off()

analysisData$predictions <- predict(mod.3, type="response")
library(reshape2)
plotData <- analysisData[,c("sample_collect_date", "location", "rolling_mean_omicron", "predictions")]
plotData <- melt(plotData,id = 1:2)
names(plotData)[3:4] <- c("Data","Concentration")

g <- guide_legend(title = "Data Source")

model.plot <- ggplot(plotData, aes(x = sample_collect_date, y = Concentration, 
                                  color = Data, shape = Data, linetype = Data)) +
  geom_point() + geom_line() + facet_wrap(~ location) + theme_bw() +
  labs(x = "Date", y = "7-day Rolling Mean Omicron Concentration (copies/µL)") + 
  scale_shape_manual(labels = c("Observed", "Model"), values = c(1,NA))  + 
  scale_linetype_manual(labels = c("Observed", "Model"), values = c("blank","solid")) + 
  scale_color_manual(labels = c("Observed", "Model"), values = c("blue","red")) +
  guides(colour = g, linetype = g, shape = g) +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle = -90, hjust=0),
        text = element_text(family = "Arial"))+
  theme_minimal()


model.plot

#save the plot
#ggsave(filename = '../Pub/Environmental Science and Technology/Figure5.jpg', plot = model.plot, width=6, height=4)