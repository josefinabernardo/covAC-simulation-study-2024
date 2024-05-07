library(tidyverse)

# Read in these files
mx_estimates <- read.csv("2024-05-01_mx_estimates_results.csv")
mx_power <- read.csv("2024-05-01_mx_power_results.csv")
gee_estimates <- read.csv("2024-05-01_gee_estimates_results.csv")
gee_power <- read.csv("2024-05-01_gee_power_results.csv")

# Appendix data
app_mx_estimates <- read.csv("2024-05-07_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-07_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-07_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-07_gee_power_appendix.csv")

# Trying to fix the data
app_gee_power$prs <- rep(c(0.02, 0.05, 0.1, 0.15), each = 9)

# Create dataset
p5_data <- app_gee_power %>%
  filter(g == 0) %>%
  select("b", "p5", 'prs') %>%
  mutate(Variable = "CT")

p6_data <- app_gee_power %>%
  filter(b == 0) %>%
  select("g", "p6", 'prs') %>%
  mutate(Variable = "SI")

colnames(p5_data)[1:3] <- c("Confounder", "Power", "PGS")
colnames(p6_data)[1:3] <- c("Confounder", "Power", "PGS")
gee_plot_data <- rbind(p5_data, p6_data)

ggplot(data = gee_plot_data, mapping = aes(x = Confounder, y = Power)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS)

