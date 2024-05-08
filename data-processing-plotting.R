library(tidyverse)

# Appendix data
app_mx_estimates <- read.csv("2024-05-08_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-08_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-08_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-08_gee_power_appendix.csv")


# Create dataset for gee results
p5_data <- app_gee_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p5", 'PGS') %>%
  mutate(Variable = "CT")

p6_data <- app_gee_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p6", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p5_data)[1:2] <- c("Confounder", "Power")
colnames(p6_data)[1:2] <- c("Confounder", "Power")
gee_plot_data <- rbind(p5_data, p6_data) %>%
  mutate(Model = "Gee")

# Create same dataset for MX
p1_data_mx <- app_mx_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p1", 'PGS') %>%
  mutate(Variable = "CT")

p2_data_mx <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p2", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p1_data_mx)[1:2] <- c("Confounder", "Power")
colnames(p2_data_mx)[1:2] <- c("Confounder", "Power")
mx_plot_data <- rbind(p1_data_mx, p2_data_mx) %>%
  mutate(Model = "OpenMx")

# Combine the two
# Merge the data frames based on Confounder, PGS, and Variable
full_plot_data <- rbind(mx_plot_data, gee_plot_data)

# Plot for only MX results
full_plot_data$PGS_percent <- scales::percent(full_plot_data$PGS)
full_plot_data$PGS_percent <- factor(full_plot_data$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

mx_plot_data$PGS_percent <- scales::percent(mx_plot_data$PGS)
mx_plot_data$PGS_percent <- factor(mx_plot_data$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

ggplot(data = mx_plot_data, mapping = aes(x = Confounder, y = Power, colour = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Variable)) +
  facet_wrap(~PGS_percent) +
  theme_light()

ggplot(data = full_plot_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  theme_light()



