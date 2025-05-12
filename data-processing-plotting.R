library(tidyverse)

# Appendix data
app_mx_estimates <- read.csv("2024-05-28_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-28_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-28_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-28_gee_power_appendix.csv")

# Debugged data
app_mx_estimates <- read.csv("debug_mx_estimates.csv")
app_mx_power <- read.csv("debug_mx_power.csv")
app_gee_estimates <- read.csv("debug_gee_estimates.csv")
app_gee_power <- read.csv("debug_gee_power.csv")


# Extract data from tables for running text

# MZ and DZ
app_mx_power %>%
  filter(round(a^2, 3) == .4, PGS == 0.1) %>%
  select(g, b, p1, p2, p3, p4, Smz, Sdz) %>%
  arrange(as.logical(b) + as.logical(g), b) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1,
         `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4) %>%
  write.table(file = pipe("pbcopy"), row.names = FALSE, quote = FALSE, sep = ",")
  #write.table(file = "clipboard", row.names = FALSE, quote = FALSE, sep = ",")

app_mx_estimates %>%
  filter(round(a^2, 3) == .4, PGS == 0.1) %>%
  select(g, b, e1, e2, e3, e4, Smz, Sdz) %>%
  arrange(as.logical(b) + as.logical(g), b) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e1,
    `SI (m2)` = e2, `CT (m3)` = e3, `SI (m3)` = e4) %>%
  write.table(file = pipe("pbcopy"), row.names = FALSE, quote = FALSE, sep = ",")
  #write.table(file = "clipboard", row.names = FALSE, quote = FALSE, sep = ",")

# DZ only
app_mx_power %>%
  filter(round(a^2, 3) == .4, PGS == 0.1) %>%
  select(g, b, p13, p14, p15, p16, Smz, Sdz) %>%
  arrange(as.logical(b) + as.logical(g), b) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p13,
         `SI (m2)` = p14, `CT (m3)` = p15, `SI (m3)` = p16) %>%
  write.table(file = pipe("pbcopy"), row.names = FALSE, quote = FALSE, sep = ",")
  #write.table(file = "clipboard", row.names = FALSE, quote = FALSE, sep = ",")

app_mx_estimates %>%
  filter(round(a^2, 3) == .4, PGS == 0.1) %>%
  select(g, b, e13, e14, e15, e16, Smz, Sdz) %>%
  arrange(as.logical(b) + as.logical(g), b) %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz)) %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e13,
         `SI (m2)` = e14, `CT (m3)` = e15, `SI (m3)` = e16) %>%
  write.table(file = pipe("pbcopy"), row.names = FALSE, quote = FALSE, sep = ",")
  #write.table(file = "clipboard", row.names = FALSE, quote = FALSE, sep = ",")

# Re-do table 4 appendix
app_mx_estimates %>%
  mutate(Smz = scales::percent(Smz), Sdz = scales::percent(Sdz),
         A = round(a^2, 1), C = round(c^2, 1), E = round(e^2, 1),
         across(e1:e9, round, 2)) %>%
  select(A, C, E, g, b, PGS, e1, e2, e3, e4, e5, e6, e7, e8, e9, Smz, Sdz) %>%
  unique() %>%
  write.table(file = pipe("pbcopy"), row.names = FALSE, quote = FALSE, sep = ",")
  #write.table(file = "clipboard", row.names = FALSE, quote = FALSE, sep = ",")


# START PLOTS

# DZ-only data
p13_mx_data <- app_mx_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  dplyr::select("b", "p13", "PGS") %>%
  mutate(Variable = "CT")

p14_mx_data <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  dplyr::select("g", "p14", 'PGS') %>%
  mutate(Variable = "SI")

colnames(p13_mx_data)[1:2] <- c("Confounder", "Power")
colnames(p14_mx_data)[1:2] <- c("Confounder", "Power")

mx_dz_data <- rbind(p13_mx_data, p14_mx_data) %>%
  mutate(Model = "OpenMx", PGS_percent = scales::percent(PGS))

# Create same dataset for Gee
p2_gee_data <- app_gee_power %>%
  filter(g == 0, round(a^2, 3) == .4) %>%
  dplyr::select("b", "p2", 'PGS') %>%
  mutate(Variable = "CT")

p3_gee_data <- app_gee_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  dplyr::select("g", "p3", 'PGS') %>%
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
  dplyr::select("b", "p1", 'PGS') %>%
  mutate(Variable = "CT")

p2_mx_data <- app_mx_power %>%
  filter(b == 0, round(a^2, 3) == .4) %>%
  dplyr::select("g", "p2", 'PGS') %>%
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


# ESTIMATES




# DETAILED PLOTS
ext_estimates <- read.csv("2024-06-04_mx_estimates_ext.csv")
ext_power <- read.csv("2024-06-04_mx_power_ext.csv")

ext_power <- ext_power %>%
                rename(`CT` = g, `SI` = b, `CT (m1)` = p1,
                        `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4)

data_noCT <- ext_power[ext_power$CT == 0, ]

# Gather the data for plotting
data_noCT_long <- gather(data_noCT, key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)

# Create the plot
ggplot(data_noCT_long, aes(x = SI, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Sibling Interaction",
       y = "Power",
       color = "Parameter") +
  theme_minimal()

# Filter the data for b = 0
data_noSI <- ext_power[ext_power$SI== 0, ]

# Gather the data for plotting
data_noSI_long <- gather(data_noSI, key = "variable", value = "value", `CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)

# Create the plot
ggplot(data_noSI_long, aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Cultural Transmission",
       y = "Power",
       color = "Parameter") +
  theme_minimal()

# Filter the data for g = 0
data_g0 <- ext_power[ext_power$g == 0, ]
data_g0_long <- gather(data_g0, key = "variable", value = "value", ,`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_g0_long$x_value <- data_g0_long$b

# Filter the data for b = 0
data_b0 <- ext_power[ext_power$b == 0, ]
data_b0_long <- gather(data_b0, key = "variable", value = "value", ,`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_b0_long$x_value <- data_b0_long$g

