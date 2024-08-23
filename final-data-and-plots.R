# Install all packages
# detach("package:gnomesims", unload = TRUE)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)
library(OpenMx)

# Run function for detailed lots in running text
paper_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                             nloci = 100,
                                             npgsloci = c(2, 5, 10, 15))

gee_data <- gnomesims::gnome_gee_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                           nloci = 100,
                                           npgsloci = c(2, 5, 10, 15))

# Create seperate data sets for processing
paper_power <- paper_data$power
paper_estimates <- paper_data$params

gee_power <- gee_data$power
gee_estimates <- gee_data$params

# Write to .csv files
write.csv(paper_estimates, file = "paper_mx_estimates.csv", row.names = TRUE)
write.csv(paper_power, file = "paper_mx_power.csv", row.names = TRUE)

write.csv(gee_estimates, file = "paper_gee_estimates.csv", row.names = TRUE)
write.csv(gee_power, file = "paper_gee_power.csv", row.names = TRUE)

# Load in the data sets
#paper_power <- read.csv("2024_08_14_mx_power.csv")
#paper_estimates <- read.csv("2024_08_14_mx_estimates.csv")

paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")

# Plot 1 - MZ & DZ Power
library(tidyverse)

p1_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p1", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p1)

p2_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p2", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Power = p2)

mx_mzdz_data <- rbind(p1_mx_data, p2_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))

# Plot 2 - MZ & DZ vs. DZ-only Power
p5_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p5", "PGS") %>%
  mutate(Variable = "CT", Sample = "DZ") %>%
  rename(Confounder = b, Power = p5)

p6_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p6", 'PGS') %>%
  mutate(Variable = "SI", Sample = "DZ") %>%
  rename(Confounder = g, Power = p6)

mx_dz_data <- rbind(p5_mx_data, p6_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

full_mx_data <- rbind(mx_mzdz_data, mx_dz_data)

ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))

# Plot 3
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
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")

# Plot 4
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
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")

# Plot 5 - old version
gee_power <- read.csv("paper_gee_power.csv")
gee_est <- read.csv("paper_gee_estimates.csv")

p1_gee_data <- gee_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p1", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p1)

p2_gee_data <- gee_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p2", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Power = p2)

gee_mzdz_data <- rbind(p1_gee_data, p2_gee_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  mutate(Method = "Gee")


mx_mzdz_data <- mx_mzdz_data %>%
  mutate(Method = "OpenMx")

full_method_data <- rbind(mx_mzdz_data, gee_mzdz_data)

ggplot(data = full_method_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))

# Debugging script

# Correlations
cor_vec <- diag(cor(paper_power[,12:19], gee_power[,12:19]))
print(round(cor_vec, 4))

# Scatterplots
vars <- paste0("p", 1:8)

for (var in vars) {
  p <- ggplot() +
    geom_point(aes_string(x = paste0("paper_power$", var), y = paste0("gee_power$", var), color = "paper_power$g")) +
    labs(
      title = paste("Scatterplot of", var, "vs", var),
      x = paste("OpenMx", var),
      y = paste("Gee", var)
    ) +
    theme_minimal() +
    scale_color_viridis_c(name = "Group")

  print(p)
}

# Scatterlots only with p1
for (var in vars) {
  p <- ggplot() +
    geom_point(aes_string(x = paste0("paper_power$", var), y = gee_power$p1, color = "paper_power$g")) +
    labs(
      title = paste("Scatterplot of", var, "vs p1"),
      x = paste("OpenMx", var),
      y = paste("Gee 1")
    ) +
    theme_minimal() +
    scale_color_viridis_c(name = "Group")

  print(p)
}


# PLot 5 - new version
p3_gee_data <- gee_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p3", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p3)

p4_gee_data <- gee_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p4", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Power = p4)

gee_combi_data <- rbind(p3_gee_data, p4_gee_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  mutate(Method = "Gee")

p3_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p3", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p3)

p4_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p4", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Power = p4)

mx_combi_data <- rbind(p3_mx_data, p4_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  mutate(Method = "OpenMx")

full_combi_data <- rbind(mx_combi_data, gee_combi_data)

ggplot(data = full_combi_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "Variable")
