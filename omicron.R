#If in Dynamics of Omicron project:
setwd("Data/") #Solana

#install.packages("tidyverse")
#install.packages("read_xcl")
#install.packages("scales")
#install.packages("plyr")
#install.packages("zoo")
#install.packages("slider")
#install.packages("ggpattern")

#Load libraries/packages
library(data.table)
library(tidyverse)
library(scales)
library(plyr)
library(dplyr)
library(zoo) # moving averages
library(slider)
library(readxl)
library(lubridate)
library(ggpubr)
library(ggpattern)

#IMPORT ALL DATA AND BEGIN SORTING

#read vaccine data and specify information about df
all.vaccines = read_excel("Latah County Vaccinations.xlsx")
colnames(all.vaccines)[3] ="sample_collect_date"
colnames(all.vaccines)[1] = 'location'
all.vaccines$sample_collect_date = as.Date(all.vaccines$sample_collect_date)
class(all.vaccines$location) = 'character'

#import case data and save dates and locations of interest
cases = read_excel("case_number_by_county.xlsx")
all.date.cases = cases %>%
  subset(select=c(sample_collect_date, Moscow_83843, Deary_83823, Genesee_83832, Potlatch_83855, Kendrick_83537, Juliaetta_83535, Troy_83871))
cases = cases %>%
  subset(sample_collect_date <= "2022-02-16" & sample_collect_date >= "2021-10-22") %>%
  subset(select=c(sample_collect_date, Moscow_83843, Deary_83823, Genesee_83832, Potlatch_83855, Kendrick_83537, Juliaetta_83535, Troy_83871))
#import omicron data
mutations = read_excel("delta_omicron_serie.xlsx") 
mutations = mutations %>% 
  subset(task == "Sample") %>% 
  separate(col = sample_id, sep = "_", into = c("location", "sample_collect_date", "RNA_ext_rep"), remove = FALSE)
mutations$sample_collect_date = as.Date(mutations$sample_collect_date)
#select date range
mutations = mutations %>%
  subset(sample_collect_date <= "2022-02-16" & sample_collect_date >= "2021-10-22")

#import N1 data
N1 = read_csv("IMCI_time_serie.csv") 
N1 <- N1 %>%
  separate(col = label, sep = "_", into = c("location", "sample_collect_date", "RNA_ext_rep", "test"), remove = FALSE)
N1 <- N1 %>%  
  subset(task == "Sample")
N1 <- N1 %>%
  subset(is.na(N1$test))
  #unite('sample_id', location:sample_collect_date, remove = FALSE) 
#extract date from sample name
N1$sample_collect_date = as.Date(N1$sample_collect_date)
#select N1 and date range, exclude other samples
N1 = N1 %>% 
  subset(pcr_gene_target == "N1") %>%
  subset(sample_collect_date <= "2022-02-16" & sample_collect_date >= "2021-10-22") %>%
  subset(!str_detect(label, "Neg")) %>%
  subset(!str_detect(label, "BCoV")) %>%
  subset(!str_detect(location, "PostFalls")) %>%
  subset(!str_detect(location, "TwinFalls")) %>%
  subset(!str_detect(location, "Deary")) %>%
  subset(!str_detect(location, "IdahoFalls")) %>%
  subset(!str_detect(location, "Rathdrum"))

N1$Variant = 'N1'

#change concentration of omicron and delta to copies/mL (was in copies/ uL)
#multiply concentration by 1000
N1$concentration <- N1$concentration * 1000

#PREPARE DELTA AND OMICRON DATA
#select delta mutations from mutations data set
delta = mutations %>%
  subset(pcr_gene_target == "L452R" | pcr_gene_target == "T478K")

#select omicron mutations fom mutations data set
omicron = mutations %>%
  subset(pcr_gene_target == "Q954H" | pcr_gene_target == "N679K")

#change concentration of omicron and delta to copies/mL (was in copies/ uL)
#multiply concentration by 1000
omicron$concentration <- omicron$concentration * 1000
delta$concentration <- delta$concentration * 1000

#calculate the mean of delta mutations for the same sample
delta.means = ddply(delta, c("sample_id", "location", "sample_collect_date"), summarise,
                    ave_concentration = mean(concentration),
                    sd = sd(concentration, na.rm = TRUE))
delta.means$Variant = 'Delta'

#calculate the mean of omicron mutations for the same sample
omicron.means = ddply(omicron, c("sample_id", "location", "sample_collect_date"), summarise,
                    ave_concentration = mean(concentration),
                    sd = sd(concentration, na.rm = TRUE))
omicron.means$Variant = 'Omicron'

#label undetected samples
omicron.means$detected = ifelse (omicron.means$ave_concentration == 0, "no", "yes")
delta.means$detected = ifelse (delta.means$ave_concentration == 0, "no", "yes")
#change concentration of omicron and delta to copies/mL (was in copies/ uL)
#multiply concentration by 1000
# Find the earliest date where omicron is greater than 0 for each location
earliest_dates <- omicron.means %>%
  filter(ave_concentration > 0) %>%
  group_by(location) %>%
  dplyr::summarize(earliest_date = min(sample_collect_date, na.rm = TRUE)) %>%
  ungroup()

#save table
#write.table(earliest_dates, "omicron_earliest_detection.tsv", sep="\t", row.names=FALSE)

#PLOT DELTA AND OMICRON INFORMATION

#define theme for plots
theme = theme(panel.background = element_rect(fill = NA),
              panel.grid.major.y = element_blank(), 
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line(colour = "grey90"), 
              panel.grid.minor.x = element_line(colour = "grey80", linetype = 3),
              axis.line = element_line(colour = "black"),
              axis.text=element_text(size=10),
              axis.text.x=element_text(angle = -90, hjust = 0),
              legend.position="top",
              legend.title=element_blank(),
              element_text(family = "Arial"))

#Plot omicron mutation means 
plot.omicron = ggplot(omicron.means, aes(x=sample_collect_date, y = ave_concentration, colour = Variant)) + 
  #ggplot(delta.means, aes(x=sample_collect_date, y = ave_concentration))
  #geom_errorbar(aes(ymax = concentration + CI, ymin = concentration - CI), width=0.15) +
  geom_point(aes(shape=detected)) + 
  scale_shape_manual(values=c(21,16))+
  geom_line() +
  facet_wrap(.~location, scales = "free") +
  theme

plot.variants = plot.omicron +
  geom_point(data = delta.means, aes(shape = detected)) + #aes(shape = as.factor(sample_type)) 
  scale_shape_manual(values=c(21,16)) +
  geom_line(data = delta.means)

#Add N1 on a secondary axis
plot.N1 = plot.variants + 
  geom_line(data = N1, aes(x=sample_collect_date, y = concentration/1), linewidth = 0.5, linetype = 2) +
  #geom_point(data = N1.N2.means, aes(x=sample_collect_date, y = ave_concentration/1, shape=detected)) +
  scale_shape_manual(values=c(21,16)) +
  scale_y_continuous(sec.axis = sec_axis(~ .*1), name = "Concentration (copies/uL in dPCR reaction)") +
  xlab(label="Collection Date") 

#save plot
#ggsave(file=paste("../Figures/omicron_delta_graphs.pdf"), plot=plot.N1, width=4, height=2.5)

#CREATE ROLLING MEAN FOR N1 AND OMICRON AND DELTA DATA SETS

#arrange data by increasing date
order.omicron.means = omicron.means[order(as.Date(omicron.means$sample_collect_date)),]


# Create a complete date sequence for each location
#complete_order.omicron <- order.omicron.means %>%
#  group_by(location) %>%
#  dplyr::summarize(min_date = min(sample_collect_date), 
#                   max_date = max(sample_collect_date), 
#                   .groups = 'drop') %>%
#  rowwise() %>%
#  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
#  select(-min_date, -max_date) %>%
#  unnest(all_dates) %>%
#  dplyr::rename(sample_collect_date = all_dates)

#  Left join the complete date sequence with the original data
#order_filled.omicron <- complete_order.omicron %>%
#  left_join(order.omicron.means, by = c("location", "sample_collect_date"))

#calculate 7-day rolling mean
#Calculate rolling averages of omicron concentrations
rolling.omicron.means = order.omicron.means %>% 
  dplyr::group_by(location) %>%
  dplyr::mutate(roll.mean = slide_index_dbl(.x = ave_concentration, 
                                     .i = sample_collect_date, 
                                     .f = ~mean(.x, na.rm = T), 
                                     .before = 6,
                                     #.after = lubridate::days(3),
                                     .complete = TRUE))
#take the average of replicates of N1 concentration
avg.N1 <- N1 %>%
  group_by(location, sample_collect_date) %>%
  dplyr::summarise(avg_concentration = mean(concentration, na.rm = TRUE))

#Calculate rolling averages of N1 concentrations
order.N1 <- avg.N1[order(as.Date(avg.N1$sample_collect_date)),]

# Create a complete date sequence for each location
complete_order.N1 <- order.N1 %>%
  group_by(location) %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
            max_date = max(sample_collect_date), 
            .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)

#  Left join the complete date sequence with the original data
order_filled.N1 <- complete_order.N1 %>%
  left_join(order.N1, by = c("location", "sample_collect_date"))
#perform 7-day rolling mean calculation for filled in all dates
rolling.N1 <- order_filled.N1 %>% 
  dplyr::group_by(location) %>%
  dplyr::mutate(roll.mean = slide_index_dbl(.x = avg_concentration, 
                                            .i = sample_collect_date, 
                                            .f = ~mean(.x, na.rm = T), 
                                            .before = 6,
                                            #.after = lubridate::days(3),
                                            .complete = TRUE))
colnames(rolling.N1)[4] <- "rolling_mean_N"
rolling.N1$data = "wastewater concentration"
#rolling.N1 <- rolling.N1[c("location", "sample_collect_date", "rolling_mean_N", "data", "variant")]

#perform 7-day rolling mean calculation for unfilled data
rolling.unfilled.N1 <- order.N1 %>% 
  dplyr::group_by(location) %>%
  dplyr::mutate(roll.mean = slide_index_dbl(.x = avg_concentration, 
                                            .i = sample_collect_date, 
                                            .f = ~mean(.x, na.rm = T), 
                                            .before = 6,
                                            #.after = lubridate::days(3),
                                            .complete = TRUE))
colnames(rolling.unfilled.N1)[4] <- "rolling_mean_N"
rolling.unfilled.N1$data = "wastewater concentration"

#Calculate rolling averages of delta concentrations
order.delta.means = delta.means[order(as.Date(delta.means$sample_collect_date)),]

# Create a complete date sequence for each location
#complete_order.delta <- order.delta.means %>%
#  group_by(location) %>%
#  dplyr::summarize(min_date = min(sample_collect_date), 
#                   max_date = max(sample_collect_date), 
#                   .groups = 'drop') %>%
#  rowwise() %>%
#  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
#  select(-min_date, -max_date) %>%
#  unnest(all_dates) %>%
#  dplyr::rename(sample_collect_date = all_dates)

#  Left join the complete date sequence with the original data
#order_filled.delta <- complete_order.delta %>%
#  left_join(order.delta.means, by = c("location", "sample_collect_date"))

#calculate 7-day rolling mean
rolling.delta.means = order.delta.means %>% 
  dplyr::group_by(location) %>%
  dplyr::mutate(roll.mean = slide_index_dbl(.x = ave_concentration, 
                                            .i = sample_collect_date, 
                                            .f = ~mean(.x, na.rm = T), 
                                            .before = 6,
                                            #.after = lubridate::days(3),
                                            .complete = TRUE))
#use pearson correlations calculated later in code to find date at which omicron is better correlated with data after 
#and delta better correlated with data before
location <- c("Moscow", "Troy", "Kendrick", "Juliaetta", "Genesee", "Potlatch")
date <- as.Date(c("2022-01-05", "2022-01-12", "2022-01-04", "2021-12-25", "2022-01-25", "2022-01-06"))
separating.date <- data.frame(location, date)

#Plot the rolling averages of the concentration of N1/N2 and omicron and delta
plot.rolling = ggplot() + 
  geom_line(data = rolling.omicron.means, aes(x=sample_collect_date, y = roll.mean, linetype = Variant, colour = Variant)) +
  #geom_point(data=omicron.means, aes(y=ave_concentration))+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week", date_labels = "%b %Y") +
  facet_wrap(.~location, scales = "free") +
  #ylim(0,3.5)+
  theme

plot.rolling.variants = plot.rolling +
  geom_line(data = rolling.delta.means, aes(x = sample_collect_date, y = roll.mean, linetype = Variant, colour = Variant))+
  #geom_point(data = delta.means, aes(y= ave_concentration))+
  xlab(label="Collection Date") +
  ylab(label="7-day Rolling Mean VOC Concentration \n (gc/ml)")
  #geom_vline(data = separating.date, mapping = aes(xintercept = date))

plot.rolling.variants
#ggsave(file=paste("../Figures/noN1_rolling_delta_omicron.pdf"), plot=plot.rolling.variants, width=20, height=12)
rolling.unfilled.N1 <- cbind(rolling.unfilled.N1, Variant = 'N1')
#Add N1 on a secondary axis
plot.rolling.all <- plot.rolling.variants + 
  geom_line(data = rolling.unfilled.N1, aes(x = sample_collect_date, y = rolling_mean_N/2, linetype = Variant, colour = Variant), size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ .*2, name = " 7-day Rolling Mean N1 Concentration (gc/ml)"))+
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
plot.rolling.all

#ggsave(file=paste("../Figures/AllVariants.jpg"), plot=plot.rolling.all, width=9, height=6, dpi = 300, units = "in")

#create shaded area showing delta and omicron peaks in different times for each location
# Define positions of vertical lines for each location
s1_vertical_lines <- data.frame(
  location = c("Moscow", "Juliaetta", "Troy","Genesee", "Potlatch"),  # Adjust these to match your locations
  xmin = ymd(c("2021-11-25", "2021-11-10", "2021-11-15","2021-12-20", "2021-11-15")),              # X-coordinate for first vertical line for spike 1
  xmax = ymd(c("2021-12-20", "2021-12-10", "2021-12-10","2022-01-20","2021-12-15"))           # X-coordinate for second vertical line for spike 1
)
s2_vertical_lines <- data.frame(
  location = c("Moscow", "Juliaetta", "Troy", "Kendrick", "Genesee", "Potlatch"),  # Adjust these to match your locations
  xmin = ymd(c("2021-12-20","2021-12-10","2022-01-01","2021-12-15","2022-01-20","2021-12-15")),      # X-coordinate for first vertical line for spike 2                     # X-coordinate for second vertical line for spike 2
  xmax = ymd(c("2022-02-16","2022-02-16","2022-02-16","2022-02-16","2022-02-16","2022-02-16"))
)
#add shaded region onto time series (can see in black-and-white but hard to look at)
plot.shaded <- plot.rolling.all +
  geom_rect_pattern(data = s1_vertical_lines, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "black", alpha = 0, pattern = 'stripe', pattern_alpha = 0.2)+
  geom_rect_pattern(data = s2_vertical_lines, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "black", alpha = 0,pattern = 'weave', pattern_alpha = 0.1)
plot.shaded
#ggsave(file=paste("../Figures/Figure2black.pdf"), plot=plot.shaded, width=7, height=5, dpi = 300, units = "in")
plot.color.blocked <- plot.rolling.all +
  geom_rect(data = s1_vertical_lines, aes(xmin = xmin, xmax = xmax, fill = "Delta"), ymin = 0, ymax = Inf, alpha = 0.4)+
  geom_rect(data = s2_vertical_lines, aes(xmin = xmin, xmax = xmax, fill ="Omicron"), ymin =0, ymax = Inf, alpha = 0.5)+
  scale_fill_manual(values=c("#fdc086","#beaed4"), labels=c("Delta", "Omicron"))+
  labs(fill = "Variant Time Period", Variant = "Variant")+
  theme(
    element_text(family = "Arial"), 
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8))

plot.color.blocked
#ggsave(file=paste("../Figures/Figure2.jpg"), plot=plot.color.blocked, width=190, height=150, dpi = 300, units = "mm")

#PREPARE CASE DATA

#remove zipcode from column names
colnames(cases) <- c('sample_collect_date','Moscow','Deary','Genesee', 'Potlatch', 'Kendrick', 'Juliaetta', 'Troy')

cases$sample_collect_date = as.Date(cases$sample_collect_date)
cases <- cbind(cases, data = "clinical cases")

#separate cases per town
Moscow.cases = cases[c("sample_collect_date", "Moscow", "data")]
colnames(Moscow.cases)[2] <- "daily_cases"

Potlatch.cases = cases[c("sample_collect_date", "Potlatch", "data")]
colnames(Potlatch.cases)[2] <- "daily_cases"

Troy.cases = cases[c("sample_collect_date", "Troy", "data")]
colnames(Troy.cases)[2] <- "daily_cases"

Genesee.cases = cases[c("sample_collect_date", "Genesee", "data")]
colnames(Genesee.cases)[2] <- "daily_cases"

Kendrick.cases = cases[c("sample_collect_date", "Kendrick", "data")]
colnames(Kendrick.cases)[2] <- "daily_cases"

Juliaetta.cases = cases[c("sample_collect_date", "Juliaetta", "data")]
colnames(Juliaetta.cases)[2] <- "daily_cases"

#create rolling means for cases for each town
order.Moscow.cases = Moscow.cases[order(as.Date(Moscow.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.M.cases <- order.Moscow.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.M.cases <- complete_order.M.cases %>%
  left_join(order.Moscow.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Moscow.cases = order_filled.M.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                            .i = sample_collect_date, 
                                            .f = ~mean(.x, na.rm = T), 
                                            .before = 6,
                                            #.after = lubridate::days(3),
                                            .complete = TRUE))

order.Troy.cases = Troy.cases[order(as.Date(Troy.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.T.cases <- order.Troy.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.T.cases <- complete_order.T.cases %>%
  left_join(order.Troy.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Troy.cases = order_filled.T.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                                     .i = sample_collect_date, 
                                                     .f = ~mean(.x, na.rm = T), 
                                                     .before = 6,
                                                     #.after = lubridate::days(3),
                                                     .complete = TRUE))

order.Juliaetta.cases = Juliaetta.cases[order(as.Date(Juliaetta.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.J.cases <- order.Juliaetta.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.J.cases <- complete_order.J.cases %>%
  left_join(order.Juliaetta.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Juliaetta.cases = order_filled.J.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                                     .i = sample_collect_date, 
                                                     .f = ~mean(.x, na.rm = T), 
                                                     .before = 6,
                                                     #.after = lubridate::days(3),
                                                     .complete = TRUE))

order.Genesee.cases = Genesee.cases[order(as.Date(Genesee.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.G.cases <- order.Genesee.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.G.cases <- complete_order.G.cases %>%
  left_join(order.Genesee.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Genesee.cases = order_filled.G.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                                     .i = sample_collect_date, 
                                                     .f = ~mean(.x, na.rm = T), 
                                                     .before = 6,
                                                     #.after = lubridate::days(3),
                                                     .complete = TRUE))

order.Kendrick.cases = Kendrick.cases[order(as.Date(Kendrick.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.K.cases <- order.Kendrick.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.K.cases <- complete_order.K.cases %>%
  left_join(order.Kendrick.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Kendrick.cases = order_filled.K.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                                     .i = sample_collect_date, 
                                                     .f = ~mean(.x, na.rm = T), 
                                                     .before = 6,
                                                     #.after = lubridate::days(3),
                                                     .complete = TRUE))

order.Potlatch.cases = Potlatch.cases[order(as.Date(Potlatch.cases$sample_collect_date)),]
# Create a complete date sequence for each location
complete_order.P.cases <- order.Potlatch.cases %>%
  dplyr::summarize(min_date = min(sample_collect_date), 
                   max_date = max(sample_collect_date), 
                   .groups = 'drop') %>%
  rowwise() %>%
  mutate(all_dates = list(seq.Date(min_date, max_date, by = "day"))) %>%
  select(-min_date, -max_date) %>%
  unnest(all_dates) %>%
  dplyr::rename(sample_collect_date = all_dates)
#  Left join the complete date sequence with the original data
order_filled.P.cases <- complete_order.P.cases %>%
  left_join(order.Potlatch.cases, by = c("sample_collect_date"))
#calculate 7-day rolling mean
rolling.Potlatch.cases = order_filled.P.cases %>% 
  dplyr::mutate(rolling_mean_cases = slide_index_dbl(.x = daily_cases, 
                                                     .i = sample_collect_date, 
                                                     .f = ~mean(.x, na.rm = T), 
                                                     .before = 6,
                                                     #.after = lubridate::days(3),
                                                     .complete = TRUE))
#create cases file with all towns 
#insert column for each town with town name
rolling.Kendrick.cases <- cbind(rolling.Kendrick.cases, location = 'Kendrick')
rolling.Juliaetta.cases <- cbind(rolling.Juliaetta.cases, location = 'Juliaetta')
rolling.Moscow.cases <- cbind(rolling.Moscow.cases, location = 'Moscow')
rolling.Troy.cases <- cbind(rolling.Troy.cases, location = 'Troy')
rolling.Potlatch.cases <- cbind(rolling.Potlatch.cases, location = 'Potlatch')
rolling.Genesee.cases <- cbind(rolling.Genesee.cases, location = 'Genesee')
#combine all df
all.towns.cases <- rbind(rolling.Kendrick.cases, rolling.Juliaetta.cases, rolling.Troy.cases, rolling.Moscow.cases, rolling.Potlatch.cases, rolling.Genesee.cases)
all.cases <- all.towns.cases %>%
  subset(sample_collect_date <= "2022-02-16" & sample_collect_date >= "2021-10-22")

#create graph with all towns split and show wastewater concentration and cases
cases.plot <- ggplot(all.cases, (aes(x=sample_collect_date, y=rolling_mean_cases, colour = data, linetype = data)))+
  geom_line()+
  labs(x = "Collection Date", y= "Clinical Cases")+
  facet_wrap(~ location, scales = "free") 
cases.plot

cases.N.plot <- cases.plot + 
  geom_line(data = rolling.N1, aes(y = rolling_mean_N/2, linetype = Variant, colour = Variant), size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ .*2, name = "N1 Concentration (copies/uL)"))+
  theme_bw(base_size = 10, ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title=element_blank())
#cases.N.plot
#ggsave(file=paste("../Figures/All_towns_N_cases.jpg"), plot=cases.N.plot, width=9.5, height=4.75, dpi = 300, units = "in")

#separate rolling N data by town 
juliaetta.N <- rolling.N1 %>%
  subset(location == "Juliaetta")

troy.N <- rolling.N1 %>%
  subset(location == "Troy")

kendrick.N <- rolling.N1 %>%
  subset(location == "Kendrick")

genesee.N <- rolling.N1 %>%
  subset(location == "Genesee")

potlatch.N <- rolling.N1 %>%
  subset(location == "Potlatch")

moscow.N = rolling.N1 %>%
  subset(location == "Moscow")

#plot rolling cases for Moscow
plot.Moscow.cases = ggplot(rolling.Moscow.cases) + 
  labs(y = NULL, x = NULL, title = 'Moscow') +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme
plot.Moscow.cases

#add rolling N values to plot with cases
plot.Moscow.cases.N1 = plot.Moscow.cases + 
  geom_line(data = moscow.N, aes(x=sample_collect_date, y = rolling_mean_N/66), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*(66)))+
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5),
        element_text(family = "Arial"), panel.grid = element_blank())

plot.Moscow.cases.N1

#ggsave(file=paste("../Figures/Moscow cases and N.jpg"), plot=plot.Moscow.cases.N1, width=5, height=4)

#plot cases and N for Potlatch
plot.Potlatch.cases = ggplot(rolling.Potlatch.cases) + 
  labs(x = NULL, y = NULL, title = "Potlatch") +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5),
        element_text(family = "Arial"), panel.grid = element_blank())

plot.Potlatch.cases

#add rolling N values to plot with cases
plot.Potlatch.cases.N1 = plot.Potlatch.cases + 
  geom_line(data = potlatch.N, aes(x=sample_collect_date, y = rolling_mean_N/500), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*500,
                                         labels = label_number(accuracy = 1),
                                         breaks = (seq(0, 10000, by =2000))))
plot.Potlatch.cases.N1
#ggsave(file=paste("../Figures/Potlatch cases and N.jpg"), plot=plot.Potlatch.cases.N1, width=5, height=4)

#plot cases and N for Genesee
plot.Genesee.cases = ggplot(rolling.Genesee.cases) + 
  labs(x = NULL, y = NULL, title = "Genesee") +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5),
        element_text(family = "Arial"), panel.grid = element_blank())
plot.Genesee.cases
#add rolling N values to plot with cases
plot.Genesee.cases.N1 = plot.Genesee.cases + 
  geom_line(data = genesee.N, aes(x=sample_collect_date, y = rolling_mean_N/3333), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*3333))
plot.Genesee.cases.N1
#ggsave(file=paste("../Figures/Genesee cases and N.jpg"), plot=plot.Genesee.cases.N1, width=5, height=4)

#plot cases and N for Juliaetta
plot.Juliaetta.cases = ggplot(rolling.Juliaetta.cases) +
  labs(x = NULL, y = NULL, title = "Juliaetta") +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5), panel.grid = element_blank())
plot.Juliaetta.cases
#add rolling N values to plot with cases
plot.Juliaetta.cases.N1 = plot.Juliaetta.cases + 
  geom_line(data = juliaetta.N, aes(x=sample_collect_date, y = rolling_mean_N/1000), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*1000))
plot.Juliaetta.cases.N1
#ggsave(file=paste("../Figures/Juliaetta cases and N.jpg"), plot=plot.Juliaetta.cases.N1, width=5, height=4)

#plot cases and N for Kendrick
plot.Kendrick.cases = ggplot(rolling.Kendrick.cases) + 
  labs(x = NULL, y = NULL, title = "Kendrick") +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5), panel.grid = element_blank())
plot.Kendrick.cases
#add rolling N values to plot with cases
plot.Kendrick.cases.N1 = plot.Kendrick.cases + 
  geom_line(data = kendrick.N, aes(x=sample_collect_date, y = rolling_mean_N/3000), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*3000, 
                                         labels = label_number(accuracy = 1),
                                         breaks = seq(0,12000, by = 2000)))
plot.Kendrick.cases.N1
#ggsave(file=paste("../Figures/Kendrick cases and N.jpg"), plot=plot.Kendrick.cases.N1, width=5, height=4)

#plot cases and N for Troy
plot.Troy.cases = ggplot(rolling.Troy.cases) +
  labs(x = NULL, y = NULL, title = "Troy") +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  geom_line(aes(x=sample_collect_date, y = rolling_mean_cases), color = "#a6611a", linetype = "solid" ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= "none", 
        plot.title = element_text(color = "black", size = 10, hjust = 0.5), panel.grid = element_blank())
plot.Troy.cases
#add rolling N values to plot with cases
plot.Troy.cases.N1 = plot.Troy.cases + 
  geom_line(data = troy.N, aes(x=sample_collect_date, y = rolling_mean_N/1000), colour = "#31a354", linetype = "dashed", size = 0.5)+ 
  scale_y_continuous(sec.axis = sec_axis(~.*1000))
plot.Troy.cases.N1
#ggsave(file=paste("../Figures/Troy cases and N.jpg"), plot=plot.Troy.cases.N1, width=5, height=4)

#create legend as a separate ggplot object
#create dataframe for legend
legend_data <- data.frame(
  x = c(1,2),
  y = c(1,2),
  label = c("Cases", "Concentration of N1"),
  linetype = c("solid", "dashed"),
  color = c("#a6611a", "#31a354")
)

#create legend
legend <- ggplot(legend_data, aes(x,y, group = label))+
  geom_line(aes(color = label, linetype = label), size = 0.7) +
  scale_color_manual(values = setNames(legend_data$color, legend_data$label))+
  scale_linetype_manual(values = setNames(legend_data$linetype, legend_data$label)) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(title = element_blank()),
    linetype = guide_legend(title = element_blank())
  )
legend

#combine plots
plot.all.cases.N1.noleg <- ggarrange(plot.Genesee.cases.N1, plot.Juliaetta.cases.N1, plot.Kendrick.cases.N1, plot.Moscow.cases.N1, plot.Potlatch.cases.N1, plot.Troy.cases.N1, nrow = 2, ncol = 3, common.legend =TRUE, legend = "bottom")
plot.all.cases.N1.noleg <- annotate_figure(plot.all.cases.N1.noleg, bottom = "Collection Date", left = text_grob("Daily New Cases", color = "#a6611a", rot = 90), right = text_grob("7-day Rolling Mean Concentration of N1 (gc/ml)", color = "#31a354", rot = 270))
plot.all.cases.N1 <- ggarrange(plot.all.cases.N1.noleg, legend, ncol = 1, heights = c(10,0.5))
plot.all.cases.N1
#ggsave(file=paste("../Figures/Figure1.jpg"), plot=plot.all.cases.N1, width=190, height=150, units = "mm", dpi = 300)

#PERFORM CORRELATION TESTS ON CLINICAL CASES AND N1 WASTEWATER CONCENTRATION

#Moscow correlation
#determine normality of data
#select needed data
mini.rolling.Moscow.cases <- rolling.Moscow.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Moscow.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(moscow.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
M.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)

#lag Moscow cases by optimal days
shifted_data <- mini.rolling.Moscow.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date - lubridate::days(5)) %>%
  select(shifted_date, rolling_mean_cases)
M.shifted <- moscow.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
M.shifted <- M.shifted[complete.cases(M.shifted$rolling_mean_N, M.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = M.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Moscow_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#Potlatch correlation
#select needed data
mini.rolling.Potlatch.cases <- rolling.Potlatch.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Potlatch.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(potlatch.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
P.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)

#lag Potlatch cases by optimal days
shifted_data <- mini.rolling.Potlatch.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(3)) %>%
  select(shifted_date, rolling_mean_cases)
P.shifted <- potlatch.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
P.shifted <- P.shifted[complete.cases(P.shifted$rolling_mean_N, P.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = P.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Potlatch_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#Troy correlation
#select needed data
mini.rolling.Troy.cases <- rolling.Troy.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Troy.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(troy.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
T.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)

#lag Troy cases by optimal days
shifted_data <- mini.rolling.Troy.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date - lubridate::days(4)) %>%
  select(shifted_date, rolling_mean_cases)
T.shifted <- troy.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
T.shifted <- T.shifted[complete.cases(T.shifted$rolling_mean_N, T.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = T.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Troy_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#Kendrick correlation
#select needed data
mini.rolling.Kendrick.cases <- rolling.Kendrick.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Kendrick.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(kendrick.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
K.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)
#lag Kendrick cases by optimal days
shifted_data <- mini.rolling.Kendrick.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date - lubridate::days(4)) %>%
  select(shifted_date, rolling_mean_cases)
K.shifted <- kendrick.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
K.shifted <- K.shifted[complete.cases(K.shifted$rolling_mean_N, K.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = K.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Kendrick_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#Genesee correlation
#select needed data
mini.rolling.Genesee.cases <- rolling.Genesee.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Genesee.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(genesee.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
G.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)
#lag Genesee cases by optimal days
shifted_data <- mini.rolling.Genesee.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date - lubridate::days(8)) %>%
  select(shifted_date, rolling_mean_cases)
G.shifted <- genesee.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
G.shifted <- G.shifted[complete.cases(G.shifted$rolling_mean_N, G.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = G.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Genesee_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#Juliaetta correlation
#select needed data
mini.rolling.Juliaetta.cases <- rolling.Juliaetta.cases[c("sample_collect_date", "rolling_mean_cases")]
#compare correlations with different lag times between clinical and wastewater data
# Create a vector to store the correlation coefficients for each lag
lag_correlations <- numeric(length = 21) # To store correlations from -10 to +10 days
p_values <- numeric(length = 21)
lags <- -10:10  # The range of lags

# Loop through each lag value
for (i in seq_along(lags)) {
  lag_value <- lags[i]
  # Create a new column for 'shifted_date' by subtracting the lag_value from 'sample_collect_date'
  shifted_case_data <- data.frame()
  shifted_case_data <- mini.rolling.Juliaetta.cases %>%
    dplyr::mutate(shifted_date = sample_collect_date + lubridate::days(lag_value)) %>%
    select(shifted_date, rolling_mean_cases)
  # Join on 'shifted_date' to align 'rolling_mean_cases' with 'rolling_mean_N' by date and delete rows with NA
  lag <- data.frame()
  lag <- dplyr::left_join(juliaetta.N, shifted_case_data, by = c("sample_collect_date" = "shifted_date"))
  lag <- lag[complete.cases(lag$rolling_mean_N, lag$rolling_mean_cases),]
  # Calculate the Pearson correlation between 'rolling_mean_N' and shifted 'rolling_mean_cases'
  correlation <- list()
  correlation <- cor.test(lag$rolling_mean_N, lag$rolling_mean_cases, method = "kendall", use = "complete.obs", exact = F)
  # Store the result in the lag_correlations vector
  lag_correlations[i] <- correlation$estimate
  p_values[i] <- correlation$p.value
}
# Create a data frame to show the lag values and their corresponding correlation coefficients
J.Kendall.results <- data.frame(Lag = lags, Correlation = lag_correlations, P_Value = p_values)
#lag Juliaetta cases by optimal days
shifted_data <- mini.rolling.Juliaetta.cases %>%
  dplyr::mutate(shifted_date = sample_collect_date - lubridate::days(6)) %>%
  select(shifted_date, rolling_mean_cases)
J.shifted <- juliaetta.N %>%
  dplyr::left_join(shifted_data, by = c("sample_collect_date" = "shifted_date"))
J.shifted <- J.shifted[complete.cases(J.shifted$rolling_mean_N, J.shifted$rolling_mean_cases),]

plot.correlation = ggscatter(data = J.shifted, x = "rolling_mean_cases", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Daily Case Count", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")

plot.correlation
#ggsave(file=paste("../Figures/Supplemental/Juliaetta_Cases_wastewater_correlation.jpg"), plot=plot.correlation)

#combine lag correlation tests for all towns
# Create an empty data frame to store all correlation results from all towns
all_results <- data.frame()
#create list with name of results
town_results_list <- list(M.Kendall.results, G.Kendall.results, P.Kendall.results, K.Kendall.results, T.Kendall.results, J.Kendall.results)
# Assuming you have a list of data frames called 'town_results_list'
# We will combine them using bind_rows()
all_results <- bind_rows(
  lapply(seq_along(town_results_list), function(i) {
    town_data <- town_results_list[[i]]
    town_data$Town <- paste("Town", i)  # Add a column to identify the town
    return(town_data)
  })
)
mean_correlations <- all_results %>%
  dplyr::group_by(Lag = as.numeric(Lag)) %>%
  dplyr::summarise(Mean_Correlation = mean(Correlation, na.rm = TRUE))
# Find the lag with the highest absolute mean correlation
best_lag <- mean_correlations %>%
  filter(abs(Mean_Correlation) == max(abs(Mean_Correlation)))

print(best_lag)
#best lag considering all results is -2

#PERFORM CORRELATION TESTS ON N1 AVG AND OMICRON AVG
#separate unfilled N1 data by location
juliaetta.unfilled.N <- rolling.unfilled.N1 %>%
  subset(location == "Juliaetta")

troy.unfilled.N <- rolling.unfilled.N1 %>%
  subset(location == "Troy")

kendrick.unfilled.N <- rolling.unfilled.N1 %>%
  subset(location == "Kendrick")

genesee.unfilled.N <- rolling.unfilled.N1 %>%
  subset(location == "Genesee")

potlatch.unfilled.N <- rolling.unfilled.N1 %>%
  subset(location == "Potlatch")

moscow.unfilled.N = rolling.unfilled.N1 %>%
  subset(location == "Moscow")

#cleanup omircon dataset
colnames(rolling.omicron.means)[8] <- "rolling_mean_omicron"
rolling.omicron.means = rolling.omicron.means[c("location", "sample_collect_date", "rolling_mean_omicron")]

#separate omicron data by town
moscow.omicron <- rolling.omicron.means %>%
  subset(location == "Moscow")

juliaetta.omicron <- rolling.omicron.means %>%
  subset(location == "Juliaetta")

troy.omicron <- rolling.omicron.means %>%
  subset(location == "Troy")

kendrick.omicron <- rolling.omicron.means %>%
  subset(location == "Kendrick")

genesee.omicron <- rolling.omicron.means %>%
  subset(location == "Genesee")

potlatch.omicron <- rolling.omicron.means %>%
  subset(location == "Potlatch")

#merge data about each town
moscow.merged <- merge(moscow.unfilled.N, moscow.omicron, by = "sample_collect_date")
moscow.merged <- moscow.merged[complete.cases(moscow.merged),]
moscow.s2.merged <- moscow.merged %>%
  subset(sample_collect_date >= "2021-12-20")
moscow.s1.merged <- moscow.merged %>%
  subset(sample_collect_date <= "2021-12-20") %>%
  subset(sample_collect_date >= "2021-11-25")

juliaetta.merged <- merge(juliaetta.unfilled.N, juliaetta.omicron, by = "sample_collect_date")
juliaetta.merged <- juliaetta.merged[complete.cases(juliaetta.merged),]
juliaetta.s2.merged <- juliaetta.merged %>%
  subset(sample_collect_date >= "2021-12-10")
juliaetta.s1.merged <- juliaetta.merged %>%
  subset(sample_collect_date <= "2021-12-10") %>%
  subset(sample_collect_date >= "2021-11-10")

troy.merged <- merge(troy.unfilled.N, troy.omicron, by = "sample_collect_date")
troy.merged <- troy.merged[complete.cases(troy.merged),]
troy.s2.merged <- troy.merged %>%
  subset(sample_collect_date >= "2022-01-01")
troy.s1.merged <- troy.merged %>%
  subset(sample_collect_date <= "2021-12-10") %>%
  subset(sample_collect_date >= "2021-11-15")

kendrick.merged <- merge(kendrick.unfilled.N, kendrick.omicron, by = "sample_collect_date")
kendrick.merged <- kendrick.merged[complete.cases(kendrick.merged),]
kendrick.s2.merged <- kendrick.merged %>%
  subset(sample_collect_date >= "2021-12-15")

genesee.merged <- merge(genesee.unfilled.N, genesee.omicron, by = "sample_collect_date")
genesee.merged <- genesee.merged[complete.cases(genesee.merged),]
genesee.s2.merged <- genesee.merged %>%
  subset(sample_collect_date >= "2022-01-20")
genesee.s1.merged <- genesee.merged %>%
  subset(sample_collect_date <= "2022-01-20") %>%
  subset(sample_collect_date >= "2021-12-20")

potlatch.merged <- merge(potlatch.unfilled.N, potlatch.omicron, by = "sample_collect_date")
potlatch.merged <- potlatch.merged[complete.cases(potlatch.merged),]
potlatch.s2.merged <- potlatch.merged %>%
  subset(sample_collect_date >= "2021-12-15")
potlatch.s1.merged <- potlatch.merged %>%
  subset(sample_collect_date <= "2021-12-15") %>%
  subset(sample_collect_date >= "2021-11-15")

#run correlation tests
#correlations for spike 1
x = moscow.s1.merged$rolling_mean_N
y = moscow.s1.merged$rolling_mean_omicron

test.result.moscow.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.moscow.s1

#correlation for spike 2
x = moscow.s2.merged$rolling_mean_N
y = moscow.s2.merged$rolling_mean_omicron
test.result.moscow.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.moscow.s2
#spike1: -0.478, pvalue = 0.120
#spike after 12-20-2021: 0.686, pvalue = 6.97e-5

plot.moscow.omi = ggscatter(data = moscow.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "kendall",
                             xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.moscow.omi
#ggsave(filename = "../Figures/Supplemental/moscow_omicron_spi2_correlation.jpg", plot = last_plot())

#correlation for spike 1
x = juliaetta.s1.merged$rolling_mean_N
y = juliaetta.s1.merged$rolling_mean_omicron

test.result.juliaetta.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.juliaetta.s1

#correlation for spike 2
x = juliaetta.s2.merged$rolling_mean_N
y = juliaetta.s2.merged$rolling_mean_omicron
test.result.juliaetta.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.juliaetta.s2
#spike 1: no omicron data
#spike after 12-05-2021: 0.679, p value=1.77e-5

plot.jul.omi = ggscatter(data = juliaetta.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "kendall",
                            xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.jul.omi
#ggsave(filename = "../Figures/Supplemental/juliaetta_omicron_correlation.jpg", plot = last_plot())

#correlation for spike 1
x = troy.s1.merged$rolling_mean_N
y = troy.s1.merged$rolling_mean_omicron

test.result.troy.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.troy.s1

#correlation for spike 2
x = troy.s2.merged$rolling_mean_N
y = troy.s2.merged$rolling_mean_omicron
test.result.troy.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.troy.s2
#spike1: no omicron data
#spike after 01-01-2022: 0.695, pvalue = 3.03e-4

plot.troy.omi = ggscatter(data = troy.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.troy.omi
#ggsave(filename = "../Figures/Supplemental/troy_omicron_correlation.jpg", plot = last_plot())

#correlation for spike 2
x = kendrick.s2.merged$rolling_mean_N
y = kendrick.s2.merged$rolling_mean_omicron
test.result.kendrick.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.kendrick.s2
#spike after 12-15-2021: 0.593, pvalue = 7.02e-4

plot.kend.omi = ggscatter(data = kendrick.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.kend.omi
#ggsave(filename = "../Figures/Supplemental/kendrick_omicron_correlation.jpg", plot = last_plot())

#correlation for spike 1
x = potlatch.s1.merged$rolling_mean_N
y = potlatch.s1.merged$rolling_mean_omicron

test.result.potlatch.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.potlatch.s1

#correlation for spike 2
x = potlatch.s2.merged$rolling_mean_N
y = potlatch.s2.merged$rolling_mean_omicron
test.result.potlatch.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.potlatch.s2
#spike1: -0.535, pvalue = 0.135
#spike after 12-15-2021:0.499, pvalue = 2.59e-3

plot.pot.omi = ggscatter(data = potlatch.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.pot.omi
#ggsave(filename = "../Figures/Supplemental/potlatch_omicron_correlation.jpg", plot = last_plot())

#correlation for spike 1
x = genesee.s1.merged$rolling_mean_N
y = genesee.s1.merged$rolling_mean_omicron

test.result.genesee.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.genesee.s1

plot.gen.omi.spi1 = ggscatter(data = genesee.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.gen.omi.spi1
#ggsave(filename = "../Figures/Supplemental/genesee_omicron_spi1_correlation.jpg", plot = last_plot())

#correlation for spike 2
x = genesee.s2.merged$rolling_mean_N
y = genesee.s2.merged$rolling_mean_omicron
test.result.genesee.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.genesee.s2

#spike1: 0.736, pvalue = 2.07e-3
#spike after 12-20-2021: 0.810, pvalue = 1.07e-2


plot.gen.omi.spi2 = ggscatter(data = genesee.s2.merged, x = "rolling_mean_omicron", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Omicron Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.gen.omi.spi2
#ggsave(filename = "../Figures/Supplemental/genesee_omicron_spi2_correlation.jpg", plot = last_plot())

#PERFORM CORRELATION TESTS ON N1 AVG AND Delta AVG

#cleanup delta dataset
colnames(rolling.delta.means)[8] <- "rolling_mean_delta"
rolling.delta.means = rolling.delta.means[c("location", "sample_collect_date", "rolling_mean_delta")]

#separate omicron data by town
moscow.delta <- rolling.delta.means %>%
  subset(location == "Moscow")

juliaetta.delta <- rolling.delta.means %>%
  subset(location == "Juliaetta")

troy.delta <- rolling.delta.means %>%
  subset(location == "Troy")

kendrick.delta <- rolling.delta.means %>%
  subset(location == "Kendrick")

genesee.delta <- rolling.delta.means %>%
  subset(location == "Genesee")

potlatch.delta <- rolling.delta.means %>%
  subset(location == "Potlatch")

#merge data about each town
moscow.dmerged <- merge(moscow.unfilled.N, moscow.delta, by = "sample_collect_date")
moscow.dmerged <- moscow.dmerged[complete.cases(moscow.dmerged),]
moscow.s2.dmerged <- moscow.dmerged %>%
  subset(sample_collect_date >= "2021-12-20")
moscow.s1.dmerged <- moscow.dmerged %>%
  subset(sample_collect_date <= "2021-12-20") %>%
  subset(sample_collect_date >= "2021-11-25")

juliaetta.dmerged <- merge(juliaetta.unfilled.N, juliaetta.delta, by = "sample_collect_date")
juliaetta.dmerged <- juliaetta.dmerged[complete.cases(juliaetta.dmerged),]
juliaetta.s2.dmerged <- juliaetta.dmerged %>%
  subset(sample_collect_date >= "2021-12-10")
juliaetta.s1.dmerged <- juliaetta.dmerged %>%
  subset(sample_collect_date <= "2021-12-10") %>%
  subset(sample_collect_date >= "2021-11-10")

troy.dmerged <- merge(troy.unfilled.N, troy.delta, by = "sample_collect_date")
troy.dmerged <- troy.dmerged[complete.cases(troy.dmerged),]
troy.s2.dmerged <- troy.dmerged %>%
  subset(sample_collect_date >= "2022-01-01")
troy.s1.dmerged <- troy.dmerged %>%
  subset(sample_collect_date <= "2021-12-10") %>%
  subset(sample_collect_date >= "2021-11-15")

kendrick.dmerged <- merge(kendrick.unfilled.N, kendrick.delta, by = "sample_collect_date")
kendrick.dmerged <- kendrick.dmerged[complete.cases(kendrick.dmerged),]
kendrick.s2.dmerged <- kendrick.dmerged %>%
  subset(sample_collect_date >= "2021-12-15")

genesee.dmerged <- merge(genesee.unfilled.N, genesee.delta, by = "sample_collect_date")
genesee.dmerged <- genesee.dmerged[complete.cases(genesee.dmerged),]
genesee.s2.dmerged <- genesee.dmerged %>%
  subset(sample_collect_date >= "2022-01-20")
genesee.s1.dmerged <- genesee.dmerged %>%
  subset(sample_collect_date <= "2022-01-20") %>%
  subset(sample_collect_date >= "2021-12-20")

potlatch.dmerged <- merge(potlatch.unfilled.N, potlatch.delta, by = "sample_collect_date")
potlatch.dmerged <- potlatch.dmerged[complete.cases(potlatch.dmerged),]
potlatch.s2.dmerged <- potlatch.dmerged %>%
  subset(sample_collect_date >= "2021-12-15")
potlatch.s1.dmerged <- potlatch.dmerged %>%
  subset(sample_collect_date <= "2021-12-15") %>%
  subset(sample_collect_date >= "2021-11-15")

#run correlation tests
#correlation for spike 1
x = moscow.s1.dmerged$rolling_mean_N
y = moscow.s1.dmerged$rolling_mean_delta

plot.moscow.del = ggscatter(data = moscow.s1.dmerged, x = "rolling_mean_delta", y = "rolling_mean_N", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "kendall",
                         xlab = "7-day Rolling Mean of Delta Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.moscow.del
#ggsave(filename = "../Figures/Supplemental/moscow_delta_correlation.jpg", plot = last_plot())

test.result.moscow.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.moscow.s1

#correlation for spike 2
x = moscow.s2.dmerged$rolling_mean_N
y = moscow.s2.dmerged$rolling_mean_delta
test.result.moscow.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.moscow.s2
#for spike1: 0.771, pvalue = 4.43e-3
#for spike2: -0.336, pvalue = 6.73e-2


#correlation for spike 1
x = juliaetta.s1.dmerged$rolling_mean_N
y = juliaetta.s1.dmerged$rolling_mean_delta

plot.juliaetta.del = ggscatter(data = juliaetta.s1.dmerged, x = "rolling_mean_delta", y = "rolling_mean_N", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "kendall",
                            xlab = "7-day Rolling Mean of Delta Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.juliaetta.del
#ggsave(filename = "../Figures/Supplemental/juliaetta_delta_correlation.jpg", plot = last_plot())

test.result.juliaetta.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.juliaetta.s1

#correlation for spike 2
x = juliaetta.s2.dmerged$rolling_mean_N
y = juliaetta.s2.dmerged$rolling_mean_delta
test.result.juliaetta.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.juliaetta.s2
#for spike1: 0.643, pvalue = 0.0260
#for spike2: 0.239, pvalue = 0.141


#correlation for spike 1
x = troy.s1.dmerged$rolling_mean_N
y = troy.s1.dmerged$rolling_mean_delta

plot.troy.del = ggscatter(data = troy.s1.dmerged, x = "rolling_mean_delta", y = "rolling_mean_N", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "kendall",
                            xlab = "7-day Rolling Mean of Delta Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.troy.del
#ggsave(filename = "../Figures/Supplemental/troy_delta_correlation.jpg", plot = last_plot())

test.result.troy.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.troy.s1

#correlation for spike 2
x = troy.s2.dmerged$rolling_mean_N
y = troy.s2.dmerged$rolling_mean_delta
test.result.troy.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.troy.s2
#for spike 1: 0.619, pvalue = 0.0509
#for spike2: -0.262, pvalue = 0.206


#correlation for spike 2
x = kendrick.s2.dmerged$rolling_mean_N
y = kendrick.s2.dmerged$rolling_mean_delta

test.result.kendrick.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.kendrick.s2
#for spike2: -0.352, pvalue = 0.0503


#correlation for spike 1
x = potlatch.s1.dmerged$rolling_mean_N
y = potlatch.s1.dmerged$rolling_mean_delta

plot.potlatch.del = ggscatter(data = potlatch.s1.dmerged, x = "rolling_mean_delta", y = "rolling_mean_N", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "kendall",
                            xlab = "7-day Rolling Mean of Delta Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.potlatch.del
#ggsave(filename = "../Figures/Supplemental/potlatch_delta_correlation.jpg", plot = last_plot())

test.result.potlatch.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.potlatch.s1

#correlation for spike 2
x = potlatch.s2.dmerged$rolling_mean_N
y = potlatch.s2.dmerged$rolling_mean_delta
test.result.potlatch.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.potlatch.s2
#for spike1: 0.333, pvalue = 0.293
#for spike2: -0.310, pvalue = 0.0934


#correlation for spike 1
x = genesee.s1.dmerged$rolling_mean_N
y = genesee.s1.dmerged$rolling_mean_delta

plot.genesee.del = ggscatter(data = genesee.s1.dmerged, x = "rolling_mean_delta", y = "rolling_mean_N", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "kendall",
                            xlab = "7-day Rolling Mean of Delta Concentration (copies/µL)", ylab = "7-day Rolling Mean of N1 Concentration (copies/µL)")
plot.genesee.del
#ggsave(filename = "../Figures/Supplemental/genesee_delta_correlation.jpg", plot = last_plot())

test.result.genesee.s1 = cor.test(x, y, method = "kendall", exact = F)
test.result.genesee.s1

#correlation for spike 2
x = genesee.s2.dmerged$rolling_mean_N
y = genesee.s2.dmerged$rolling_mean_delta
test.result.genesee.s2 = cor.test(x, y, method = "kendall", exact = F)
test.result.genesee.s2

#for spike 1: 0.609, pvalue = 0.0126
#for spike 2: no delta data, all 0


#combine df so there is 1 df with delta, N1 and omicron for each town
moscow.merged <- merge(moscow.merged, moscow.delta, by = "sample_collect_date")
troy.merged <- merge(troy.merged, troy.delta, by = "sample_collect_date")
genesee.merged <- merge(genesee.merged, genesee.delta, by = "sample_collect_date")
kendrick.merged <- merge(kendrick.merged, kendrick.delta, by = "sample_collect_date")
potlatch.merged <- merge(potlatch.merged, potlatch.delta, by = "sample_collect_date")
juliaetta.merged <- merge(juliaetta.merged, juliaetta.delta, by = "sample_collect_date")

#cleanup merged df
juliaetta.merged <- select(juliaetta.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))
moscow.merged <- select(moscow.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))
kendrick.merged <- select(kendrick.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))
troy.merged <- select(troy.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))
potlatch.merged <- select(potlatch.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))
genesee.merged <- select(genesee.merged, c('sample_collect_date', 'location.x', 'rolling_mean_N', 'rolling_mean_omicron', 'rolling_mean_delta'))

#rename merged columns
colnames(juliaetta.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")
colnames(troy.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")
colnames(moscow.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")
colnames(kendrick.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")
colnames(potlatch.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")
colnames(genesee.merged) <- c("sample_collect_date", "location", "rolling_mean_N", "rolling_mean_omicron", "rolling_mean_delta")

#add column with population data
juliaetta.merged <- cbind(juliaetta.merged, "population" = 1167)
moscow.merged <- cbind(moscow.merged, "population" = 26739)
troy.merged <- cbind(troy.merged, "population" = 2015)
kendrick.merged <- cbind(kendrick.merged, "population" = 985)
potlatch.merged <- cbind(potlatch.merged, "population" = 2115)
genesee.merged <- cbind(genesee.merged, "population" = 1701)

#combine into one 
all.wastewater <- rbind(juliaetta.merged, moscow.merged, troy.merged, kendrick.merged, potlatch.merged, genesee.merged)

#save wastewater data
#write.table(all.wastewater, "wastewater_data.tsv", sep="\t", row.names=FALSE)

#PREPARE VACCINE DATA

#rename locations
all.vaccines$location <- str_replace(all.vaccines$location, '83843', 'Moscow')
all.vaccines$location <- str_replace(all.vaccines$location, '83535', 'Juliaetta')
all.vaccines$location <- str_replace(all.vaccines$location, '83537', 'Kendrick')
all.vaccines$location <- str_replace(all.vaccines$location, '83871', 'Troy')
all.vaccines$location <- str_replace(all.vaccines$location, '83855', 'Potlatch')
all.vaccines$location <- str_replace(all.vaccines$location, '83832', 'Genesee')

#select zipcodes we care about, count up how many of each dose per day and zipcode, and combine replicate rows
vaccines = all.vaccines %>% 
  subset(location == "Troy" | location == "Moscow" | location == "Juliaetta" | location == "Kendrick" | 
           location == "Genesee" | location == "Potlatch") %>%
  group_by(location, sample_collect_date, DoseNumber) %>%
  dplyr::mutate(Number_of_Vaccines = n()) %>%
  distinct() %>%
  ungroup()

#fill in missing dates for all location and DoseNumber combinations 
# with a 0 value for number of vaccines
# Create a sequence of dates from 2020-12-15 to 2023-01-04
all_dates <- seq(as.Date("2020-12-15"), as.Date("2023-01-04"), by = "day")

# Get all unique combinations of location and DoseNumber from your data
location_dose_combinations <- vaccines %>%
  select(location, DoseNumber) %>%
  distinct()

# Create a data frame of all combinations of location, DoseNumber, and all dates
expanded_grid <- location_dose_combinations %>%
  crossing(sample_collect_date = all_dates)

# Merge this expanded grid with the original data
vax.result <- expanded_grid %>%
  left_join(vaccines, by = c("location", "DoseNumber", "sample_collect_date")) %>%
  mutate(Number_of_Vaccines = if_else(is.na(Number_of_Vaccines), 0, Number_of_Vaccines))

#remove county column from df
vax.result <- select(vax.result, c("location", "DoseNumber", "sample_collect_date", "Number_of_Vaccines"))

#standardize vaccine data by percent of population
Moscow.vax <- vaccines %>%
  subset(location == 'Moscow')
Moscow.vax$StandardizedVaccines <- (Moscow.vax$Number_of_Vaccines*100)/26739

Troy.vax <- vaccines %>%
  subset(location == 'Troy')
Troy.vax$StandardizedVaccines <- (Troy.vax$Number_of_Vaccines*100)/2015

Juliaetta.vax <- vaccines %>%
  subset(location == 'Juliaetta')
Juliaetta.vax$StandardizedVaccines <- (Juliaetta.vax$Number_of_Vaccines*100)/1167

Kendrick.vax <- vaccines %>%
  subset(location == 'Kendrick')
Kendrick.vax$StandardizedVaccines <- (Kendrick.vax$Number_of_Vaccines*100)/985

Genesee.vax <- vaccines %>%
  subset(location == 'Genesee')
Genesee.vax$StandardizedVaccines <- (Genesee.vax$Number_of_Vaccines*100)/1701

Potlatch.vax <- vaccines %>%
  subset(location == 'Potlatch')
Potlatch.vax$StandardizedVaccines <- (Potlatch.vax$Number_of_Vaccines*100)/2115

standardized.vax <- rbind(Moscow.vax,Potlatch.vax,Troy.vax,Genesee.vax,Juliaetta.vax,Kendrick.vax)

#separate data to just look at individual shots
shot.1 <- standardized.vax %>%
  subset(DoseNumber == '1')
shot.1 <- cbind(shot.1,labels='1st Dose')
shot.2 <- standardized.vax %>%
  subset(DoseNumber == '2')
shot.2 <- cbind(shot.2,labels='Fully Vaccinated')
shot.3 <- standardized.vax %>%
  subset(DoseNumber == '3')
shot.3 <- cbind(shot.3,labels='1st Booster')
shot.4 <- standardized.vax %>%
  subset(DoseNumber == '4')
shot.4 <- cbind(shot.4,labels='2nd Booster')

#take the cumulative sum of the standardized vaccines per shot
standardized.vax <- rbind(shot.1,shot.2, shot.3,shot.4)

standardized.sum.vaccines = standardized.vax %>%
  group_by(location, DoseNumber) %>%
  dplyr::mutate(TotalVaccines = cumsum(StandardizedVaccines))

#create a plot showing the uptake of each vaccine, shown as a cumulative sum curve
plot.each.shot = ggplot(standardized.sum.vaccines, aes(x=sample_collect_date, y = TotalVaccines, colour = labels, shape = labels))+
  geom_point()+
  geom_vline(xintercept = as.Date("2021-12-10"), linetype = "solid", color = "black")+
  geom_vline(xintercept = as.Date("2021-06-10"), linetype = "dashed", color = "black")+
  xlab('Date')+
  ylim(0,70)+
  ylab('Percent of Population Vaccinated')+
  facet_wrap(.~location, scales = "free")+
  theme_bw(base_size = 10) +
  theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(angle = -90, hjust = 0), panel.grid = element_blank()) +
  #changes order of legend labels for both color and shape
  scale_color_discrete(limits = (c('1st Dose', 'Fully Vaccinated', '1st Booster', '2nd Booster'))) +
  scale_shape_discrete(limits = (c('1st Dose', 'Fully Vaccinated', '1st Booster', '2nd Booster'))) 
plot.each.shot
#ggsave(file=paste("../Figures/Figure3.jpg"), plot=plot.each.shot, width=190, height=150, units = "mm", dpi = 300)

#create vaccine dataframe with counts per day and columns per shot
#separate vax.result by shot to create columns shot 1, shot 2, shot 3 and shot 4 of Number_of_Vaccines
vax1.num <- vax.result %>%
  subset(DoseNumber == '1')
vax2.num <- vax.result %>%
  subset(DoseNumber == '2')
vax3.num <- vax.result %>%
  subset(DoseNumber == '3')
vax4.num <- vax.result %>%
  subset(DoseNumber == '4')

#remove unneeded columns
vax1.num <- select(vax1.num, c("location", "sample_collect_date", "Number_of_Vaccines"))
vax2.num <- select(vax2.num, c("location", "sample_collect_date", "Number_of_Vaccines"))
vax3.num <- select(vax3.num, c("location", "sample_collect_date", "Number_of_Vaccines"))
vax4.num <- select(vax4.num, c("location", "sample_collect_date", "Number_of_Vaccines")) 

#rename TotalVaccines column
colnames(vax1.num) <- c("location", "sample_collect_date", "shot_1")
colnames(vax2.num) <- c("location", "sample_collect_date", "shot_2")
colnames(vax3.num) <- c("location", "sample_collect_date", "shot_3")
colnames(vax4.num) <- c("location", "sample_collect_date", "shot_4")

#combine all vaccine sum dataframes
clean.num.vaccines <- merge(vax1.num, vax2.num, by = c("location", "sample_collect_date"), all = TRUE)
clean.num.vaccines <- merge(clean.num.vaccines, vax3.num, by = c("location", "sample_collect_date"), all = TRUE)
clean.num.vaccines <- merge(clean.num.vaccines, vax4.num, by = c("location", "sample_collect_date"), all = TRUE)

#delete unneeded columns
clean.num.vaccines <- select(clean.num.vaccines, c("location", "sample_collect_date", "shot_1", "shot_2", "shot_3", "shot_4"))

#save number of vaccines per location per shot per day
#write.table(clean.num.vaccines, "vaccines_data.tsv", sep="\t", row.names=FALSE)

#separate standardized.sum.vaccines by shot to create columns shot 1, shot 2, shot 3 and shot 4 of TotalVaccines
sum.1 <- standardized.sum.vaccines %>%
  subset(DoseNumber == '1')
sum.2 <- standardized.sum.vaccines %>%
  subset(DoseNumber == '2')
sum.3 <- standardized.sum.vaccines %>%
  subset(DoseNumber == '3')
sum.4 <- standardized.sum.vaccines %>%
  subset(DoseNumber == '4')

#remove unneeded columns
sum.1 <- select(sum.1, c("location", "sample_collect_date", "TotalVaccines"))
sum.2 <- select(sum.2, c("location", "sample_collect_date", "TotalVaccines"))
sum.3 <- select(sum.3, c("location", "sample_collect_date", "TotalVaccines"))
sum.4 <- select(sum.4, c("location", "sample_collect_date", "TotalVaccines")) 

#rename TotalVaccines column
colnames(sum.1) <- c("DoseNumber", "location", "sample_collect_date", "shot_1")
colnames(sum.2) <- c("DoseNumber", "location", "sample_collect_date", "shot_2")
colnames(sum.3) <- c("DoseNumber", "location", "sample_collect_date", "shot_3")
colnames(sum.4) <- c("DoseNumber", "location", "sample_collect_date", "shot_4")

#combine all vaccine sum dataframes
better.vaccines <- merge(sum.1, sum.2, by = c("location", "sample_collect_date"), all = TRUE)
better.vaccines <- merge(better.vaccines, sum.3, by = c("location", "sample_collect_date"), all = TRUE)
better.vaccines <- merge(better.vaccines, sum.4, by = c("location", "sample_collect_date"), all = TRUE)

#delete unneeded columns
better.vaccines <- select(better.vaccines, c("location", "sample_collect_date", "shot_1", "shot_2", "shot_3", "shot_4"))

#save vaccines df
#write.table(better.vaccines, "cumulative_vaccines_data.tsv", sep="\t", row.names=FALSE)

#delete unneeded columns in cases df
better.cases <- select(all.cases, c("location", "sample_collect_date", "rolling_mean_cases"))
all.better.cases <- select(all.towns.cases, c("location", "sample_collect_date", "rolling_mean_cases"))
#save cases dataframe
#write.table(better.cases, "cases_data.tsv", sep="\t", row.names=FALSE)
#write.table(all.better.cases, "all_cases_data.tsv", sep="\t", row.names=FALSE)

