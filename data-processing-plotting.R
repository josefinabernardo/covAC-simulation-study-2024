library(tidyverse)

# Appendix data
app_mx_estimates <- read.csv("2024-05-08_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-08_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-08_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-08_gee_power_appendix.csv")

# POWER

# Create dataset for gee results
p6_data <- app_gee_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p6", 'PGS') %>%
  mutate(Variable = "CT")

p7_data <- app_gee_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p7", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p6_data)[1:2] <- c("Confounder", "Power")
colnames(p7_data)[1:2] <- c("Confounder", "Power")
gee_plot_data <- rbind(p6_data, p7_data) %>%
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

# ESTIMATES

# Create dataset for gee results
e6_data <- app_gee_estimates %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "e6", 'PGS') %>%
  mutate(Variable = "CT")

e7_data <- app_gee_estimates %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "e7", 'PGS') %>%
  mutate(Variable = "SI")

colnames(e6_data)[1:2] <- c("Confounder", "Power")
colnames(e7_data)[1:2] <- c("Confounder", "Power")
gee_plot_data_est <- rbind(e6_data, e7_data) %>%
  mutate(Model = "Gee")

# Create same dataset for MX
e1_data_mx <- app_mx_estimates %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "e1", 'PGS') %>%
  mutate(Variable = "CT")

e2_data_mx <- app_mx_estimates %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "e2", 'PGS') %>%
  mutate(Variable = "SI")

colnames(e1_data_mx)[1:2] <- c("Confounder", "Power")
colnames(e2_data_mx)[1:2] <- c("Confounder", "Power")
mx_plot_data_est <- rbind(e1_data_mx, e2_data_mx) %>%
  mutate(Model = "OpenMx")

# Combine the two
# Merge the data frames based on Confounder, PGS, and Variable
full_plot_data_est <- rbind(mx_plot_data_est, gee_plot_data_est)

# Plot for only MX results
full_plot_data_est$PGS_percent <- scales::percent(full_plot_data_est$PGS)
full_plot_data_est$PGS_percent <- factor(full_plot_data_est$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

mx_plot_data_est$PGS_percent <- scales::percent(mx_plot_data_est$PGS)
mx_plot_data_est$PGS_percent <- factor(mx_plot_data_est$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

ggplot(data = mx_plot_data_est, mapping = aes(x = Confounder, y = Power, colour = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Variable)) +
  facet_wrap(~PGS_percent) +
  theme_light()

ggplot(data = full_plot_data_est, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  theme_light()

# Compare MZ and DZ to DZ only

dz_p13_data <- app_mx_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p13", 'PGS') %>%
  mutate(Variable = "CT")

dz_p14_data <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p14", 'PGS') %>%
  mutate(Variable = "SI")

colnames(dz_p13_data)[1:2] <- c("Confounder", "Power")
colnames(dz_p14_data)[1:2] <- c("Confounder", "Power")

dz_plot_data <- rbind(dz_p13_data, dz_p14_data) %>%
  mutate(Sample = "DZ", PGS_percent = scales::percent(PGS))

mz_plot_data <- mx_plot_data %>%
  select(-Model,) %>%
  mutate(Sample = "MZ & DZ")

mzdz_plot_data <- rbind(mz_plot_data, dz_plot_data)

ggplot(data = mzdz_plot_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  theme_light()
