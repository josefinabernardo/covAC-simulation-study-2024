library(tidyverse)

# Appendix data
app_mx_estimates <- read.csv("2024-05-08_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-08_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-08_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-08_gee_power_appendix.csv")

# DZ-only data
p13_mx_data <- app_mx_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p13", 'PGS') %>%
  mutate(Variable = "CT")

p14_mx_data <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p14", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p13_mx_data)[1:2] <- c("Confounder", "Power")
colnames(p14_mx_data)[1:2] <- c("Confounder", "Power")

mx_dz_data <- rbind(p13_mx_data, p14_mx_data) %>%
  mutate(Model = "OpenMx", PGS_percent = scales::percent(PGS))

# Create same dataset for MX
p2_gee_data <- app_gee_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p2", 'PGS') %>%
  mutate(Variable = "CT")

p3_gee_data <- app_gee_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p3", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p2_gee_data)[1:2] <- c("Confounder", "Power")
colnames(p3_gee_data)[1:2] <- c("Confounder", "Power")
gee_dz_data <- rbind(p2_gee_data, p3_gee_data) %>%
  mutate(Model = "Gee", PGS_percent = scales::percent(PGS))

# Combine the two
# Merge the data frames based on Confounder, PGS, and Variable
full_dz_data <- rbind(mx_dz_data, gee_dz_data)

full_dz_data$PGS_percent <- scales::percent(full_dz_data$PGS)
full_dz_data$PGS_percent <- factor(full_dz_data$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

full_dz_data %>%
  filter(Model == "OpenMx") %>%
  ggplot(data = ., mapping = aes(x = Confounder, y = Power, colour = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Variable)) +
  facet_wrap(~PGS_percent) +
  theme_light()

ggplot(data = full_dz_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  theme_light()

# Plot of OpenMx MZ&DZ vs. DZ-only
p1_mx_data <- app_mx_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  select("b", "p1", 'PGS') %>%
  mutate(Variable = "CT")

p2_mx_data <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  select("g", "p2", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p1_mx_data)[1:2] <- c("Confounder", "Power")
colnames(p2_mx_data)[1:2] <- c("Confounder", "Power")

mx_mzdz_data <- rbind(p1_mx_data, p2_mx_data) %>%
  mutate(Model = "OpenMx", Sample = "MZ & DZ", PGS_percent = scales::percent(PGS))

mx_dz_data <- mx_dz_data %>%
  mutate(Sample = "DZ")

# Plot OpenMx DZ vs. MZ & DZ
full_mx_data <- rbind(mx_dz_data, mx_mzdz_data)

full_mx_data$PGS_percent <- scales::percent(full_mx_data$PGS)
full_mx_data$PGS_percent <- factor(full_mx_data$PGS_percent, levels = c("2%", "5%", "10%", "15%"))

ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  theme_light()
