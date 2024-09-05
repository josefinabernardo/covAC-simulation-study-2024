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
                                           npgsloci = c(2, 5, 10, 15), cmethod = "exchangeable")

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
#paper_

# Start here if files are already created
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

plot1 <- ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")

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

plot2 <- ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")

# Plot 3
ext_estimates <- read.csv("2024-06-04_mx_estimates_ext.csv")
ext_power <- read.csv("2024-06-04_mx_power_ext.csv")

ext_power <- ext_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1,
         `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4)

data_noCT <- ext_power[ext_power$CT == 0, ]

# Gather the data for plotting
data_noCT_long <- gather(data_noCT, key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_noCT_long$variable <- factor(data_noCT_long$variable,
                                  levels = c("CT (m1)", "SI (m2)", "CT (m3)", "SI (m3)"))

# Create the plot
plot3 <- ggplot(data_noCT_long, aes(x = SI, y = value, color = variable)) +
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
data_noSI_long$variable <- factor(data_noSI_long$variable,
                                  levels = c("CT (m1)", "SI (m2)", "CT (m3)", "SI (m3)"))

# Create the plot
plot4 <- ggplot(data_noSI_long, aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Cultural Transmission",
       y = "Power",
       color = "Parameter") +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")


# Plot 5 - new version
plot5_mx_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                             nloci = 100,
                                             npgsloci = c(2, 5, 10, 15))

plot5_gee_data <- gnomesims::gnome_gee_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                            nloci = 100,
                                            npgsloci = c(2, 5, 10, 15), cmethod = "exchangeable")

plot5_mx_power <- plot5_mx_data$power
plot5_mx_estimates <- plot5_mx_data$params

plot5_gee_power <- plot5_gee_data$power
plot5_gee_estimates <- plot5_gee_data$params

write.csv(plot5_mx_power, "plot5_mx_power_estimates.csv", row.names = TRUE)
write.csv(plot5_mx_estimates, "plot5_mx_parameter_estimates.csv", row.names = TRUE)
write.csv(plot5_gee_power, "plot5_gee_power_estimates.csv", row.names = TRUE)
write.csv(plot5_gee_estimates, "plot5_gee_parameter_estimates.csv", row.names = TRUE)

mx_power5 <- read.csv("plot5_mx_power_estimates.csv")
gee_power5 <- read.csv("plot5_gee_power_estimates.csv")

p3_gee_data <- gee_power5 %>%
#  dplyr::filter(PGS == .10) %>%
  dplyr::select("g", "b", "p3", "PGS") %>%
  mutate(Variable = "CT") %>%
  rename("True Parameter" = "g", "Other covAC Source" = "b", "Power" = "p3")

p4_gee_data <- gee_power5 %>%
#  dplyr::filter(PGS == .10) %>%
  dplyr::select("b", "g", "p4", "PGS") %>%
  mutate(Variable = "SI") %>%
  rename("True Parameter" = "b", "Other covAC Source" = "g", "Power" = "p4")

gee_combi_data <- rbind(p3_gee_data, p4_gee_data) %>%
  mutate(Method = "Gee",
         Label = paste0(Method, " (Other Source: ", `Other covAC Source`, ")"))

p3_mx_data <- mx_power5 %>%
#  dplyr::filter(PGS == .10) %>%
  dplyr::select("g", "b", "p3", "PGS") %>%
  mutate(Variable = "CT") %>%
  rename("True Parameter" = "g", "Other covAC Source" = "b", "Power" = "p3")

p4_mx_data <- mx_power5 %>%
#  dplyr::filter(PGS == .10) %>%
  dplyr::select("b", "g", "p4", "PGS") %>%
  mutate(Variable = "SI") %>%
  rename("True Parameter" = "b", "Other covAC Source" = "g", "Power" = "p4")

mx_combi_data <- rbind(p3_mx_data, p4_mx_data) %>%
  mutate(Method = "OpenMx",
         Label = paste0(Method, " (Other Source: ", `Other covAC Source`, ")"))

full_combi_data <- rbind(mx_combi_data, gee_combi_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

ct_data <- full_combi_data %>%
  dplyr::filter(Variable == "CT", `Other covAC Source` %in% c(0, .05, .1))


si_data <- full_combi_data %>%
  dplyr::filter(Variable == "SI", `Other covAC Source` %in% c(0, .05, .1))

ggplot(data = ct_data, mapping = aes(x = `True Parameter`, y = Power, group = as.factor(Label), color = as.factor(`Other covAC Source`))) +
  geom_point() +
  geom_line(aes(linetype = Method)) +
  facet_wrap(~PGS_percent) +
  theme(text = element_text(family = "serif")) +
  scale_linetype_manual(values = c("OpenMx" = "solid", "Gee" = "dotted")) +
  scale_color_viridis_d(name = "True Value Sibling Interaction (b)") +
  labs(x = "True Value Cultural Transmission (g)", y = "Power") +
  theme_light() +
  theme(text = element_text(family = "serif"))

ggplot(data = si_data, mapping = aes(x = `True Parameter`, y = Power, group = as.factor(Label), color = as.factor(`Other covAC Source`))) +
  geom_point() +
  geom_line(aes(linetype = Method)) +
  facet_wrap(~PGS_percent) +
  theme(text = element_text(family = "serif")) +
  scale_linetype_manual(values = c("OpenMx" = "solid", "Gee" = "dotted")) +
  scale_color_viridis_d(name = "True Value Cultural Transmission (g)") +
  labs(x = "True Value Sibling Interaction (b)", y = "Power") +
  theme_light() +
  theme(text = element_text(family = "serif"))

# Select a single PGS Level
full_combi_data %>%
  mutate(Label = paste(Method, `Other covAC Source`),
         Variable = recode(Variable, "CT" = "Cultural Transmission (CT)", "SI" = "Sibling Interaction (SI)")) %>%
  filter(PGS == .1, `Other covAC Source` %in% c(0, .05, .1)) %>%
ggplot(data = ., mapping = aes(x = `True Parameter`, y = Power,
                               color = as.factor(`Other covAC Source`), linetype = Method,
                               group = Label)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Variable) +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "Other Source of covAC") +
  scale_linetype_manual(values = c("OpenMx" = "solid", "Gee" = "dotted")) +
  labs(x = "True Parameter (g for CT, b for SI)", y = "Power") +
  theme_light() +
  theme(text = element_text(family = "serif"))


full_combi_data %>%
  mutate(Label = paste(Method, `Other covAC Source`, PGS_percent),
         Variable = recode(Variable, "CT" = "Cultural Transmission (CT)", "SI" = "Sibling Interaction (SI)")) %>%
  filter(PGS_percent %in% c("2%", "10%", "15%"),
         `Other covAC Source` %in% c(0, .05, .1)) %>%
  ggplot(aes(x = `True Parameter`, y = Power,
             color = as.factor(`Other covAC Source`), linetype = Method,
             group = Label)) +
  geom_point() +
  geom_line() +
  facet_grid(PGS_percent ~ Variable) +  # Custom facet labels
  scale_color_viridis_d(name = "Other Source of covAC") +
  scale_linetype_manual(values = c("OpenMx" = "solid", "Gee" = "dotted"), name = "Method") +  # Custom legend title
  labs(x = "True Parameter (g for CT, b for SI)", y = "Power") +
  theme_light() +
  theme(text = element_text(family = "serif"))

# Combine CT and SI into one plot
full_combi_data %>%
  filter(`Other covAC Source` %in% c(.1)) %>%
  mutate(Label = paste(Method, Variable)) %>%
  ggplot(data = ., mapping = aes(x = `True Parameter`, y = Power,
                                 group = as.factor(Label), linetype = Method,
                                 colour = Variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(~PGS_percent) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "Variable") +
  scale_linetype_manual(values = c("OpenMx" = "solid", "Gee" = "dotted"))

# Set width and height in inches (1 inch = 25.4 mm)
width <- 200 / 25.4  # Convert 200 mm to inches
height <- 150 / 25.4  # Convert 150 mm to inches

# Export each plot as a PDF with the desired size
ggsave("plot1.pdf", plot = plot1, device = "pdf", width = width, height = height)
ggsave("plot2.pdf", plot = plot2, device = "pdf", width = width, height = height)
ggsave("plot3.pdf", plot = plot3, device = "pdf", width = width, height = height)
ggsave("plot4.pdf", plot = plot4, device = "pdf", width = width, height = height)
ggsave("plot5.pdf", plot = plot5, device = "pdf", width = width, height = height)

