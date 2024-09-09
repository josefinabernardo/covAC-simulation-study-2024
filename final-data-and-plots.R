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


# Plot 5
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
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

gee_mzdz_data$Method <- c("Gee")
mx_mzdz_data$Method <- c("OpenMx")

full_mzdz_data <- rbind(gee_mzdz_data, mx_mzdz_data)

plot5 <- ggplot(data = full_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
geom_line(linewidth = 0.8) +
geom_point() +
facet_wrap(~PGS_percent) +
scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
jtools::theme_apa() +
theme(text = element_text(family = "serif")) +
scale_color_viridis_d(name = "variable")


# Set width and height in inches (1 inch = 25.4 mm)
width <- 200 / 25.4  # Convert 200 mm to inches
height <- 150 / 25.4  # Convert 150 mm to inches

# Export each plot as a PDF with the desired size
ggsave("plot1.pdf", plot = plot1, device = "pdf", width = width, height = height)
ggsave("plot2.pdf", plot = plot2, device = "pdf", width = width, height = height)
ggsave("plot3.pdf", plot = plot3, device = "pdf", width = width, height = height)
ggsave("plot4.pdf", plot = plot4, device = "pdf", width = width, height = height)
ggsave("plot5.pdf", plot = plot5, device = "pdf", width = width, height = height)

ggplot(data = full_method_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d()

# Calculate PGS Predictive Power
gnome_pgs_pred <- function(a, c, e, g, b, pgs_prop, vA = 1, vC = 1, vE = 1) {

  # Phenotypic variances
  var_mz = 2*(g+a/2+b/2)**2*vA+(a*vA/2+b*vA/2)*a+(a*vA/2+b*vA/2)*b+c**2*vC+e**2*vE
  var_dz = 2*(g+a/2+b/2)**2*vA+a**2*vA/2+b**2*vA/2+c**2*vC+e**2*vE

  # PGS variance
  vA1 = (1-pgs_prop)*vA
  vA2  = pgs_prop*vA

  var_pgs = 2*(g+a/2+b/2)**2 *vA2 +a**2*vA2 /2+b**2*vA2 /2

  effect_mz <- var_pgs/var_mz
  effect_dz <- var_pgs/var_dz

  effect = list(mz = effect_mz, dz = effect_dz)
  return(effect)
}

# Add these variable to
paper_power %>%
  filter(b == 0, g == 0) %>%
  mutate(PGS_effect_MZ = gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$mz,
         PGS_effect_DZ = gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$dz)

paper_power %>%
  filter(PGS == 0.10) %>%
  mutate(PGS_effect_MZ = gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$mz,
         PGS_effect_DZ = gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$dz)

appendix_table5 <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.05), si = seq(0,.1,.05),
                               nloci = 100,
                               npgsloci = 10)

appendix_table5$power %>%
  mutate(PGS_effect_MZ = scales::percent(round(gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$mz, 4)),
         PGS_effect_DZ = scales::percent(round(gnome_pgs_pred(a, c, e, g, b, pgs_prop = PGS)$dz, 4))) %>%
  select(g, b, PGS_effect_MZ, PGS_effect_DZ)

appendix_table5$power %>%
  mutate(Smz = scales::percent(round(Smz,4)),
         Sdz = scales::percent(round(Sdz,4)),
         g = round(g, 3), b = round(b, 3)) %>%
  select(g, b, p1:p4, Smz, Sdz)

appendix_table5$params %>%
  mutate(Smz = scales::percent(round(Smz,4)),
         Sdz = scales::percent(round(Sdz,4)),
         g = round(g, 3), b = round(b, 3)) %>%
  select(g, b, e1:e4, Smz, Sdz)

# Final tables gee
gee_appendix_table_data <- gnomesims::gnome_gee_simulation(ct = seq(0,.1,.05), si = seq(0,.1,.05),
                                            nloci = 100,
                                            npgsloci = c(2, 5, 10, 15), cmethod = "exchangeable")

gee_appendix_table3 <- gee_appendix_table_data$power %>%
  mutate(Smz = scales::percent(Smz, accuracy = 1),
         Sdz = scales::percent(Sdz, accuracy = 1),
         g = round(g, 2),
         b = round(b, 2),
         A = round(a^2, 1),
         C = round(c^2, 1),
         E = round(e^2, 1),
         p1 = round(p1, 2),
         p2 = round(p2, 2),
         p3 = round(p3, 2),
         p4 = round(p4, 2),
         p5 = round(p5, 2),
         p6 = round(p6, 2),
         p7 = round(p7, 2),
         p8 = round(p8, 2)) %>%
  select(A, C, E, g, b, PGS, p1:p8, Smz, Sdz)


gee_appendix_table4 <- gee_appendix_table_data$params %>%
  mutate(Smz = scales::percent(Smz, accuracy = 1),
         Sdz = scales::percent(Sdz, accuracy = 1),
         g = round(g, 2),
         b = round(b, 2),
         A = round(a^2, 1),
         C = round(c^2, 1),
         E = round(e^2, 1),
         e1 = round(e1, 2),
         e2 = round(e2, 2),
         e3 = round(e3, 2),
         e4 = round(e4, 2),
         e5 = round(e5, 2),
         e6 = round(e6, 2),
         e7 = round(e7, 2),
         e8 = round(e8, 2)) %>%
  select(A, C, E, g, b, PGS, e1:e8, Smz, Sdz)

# Convert gee_appendix_table3 to CSV format without index and print
cat(paste0(capture.output(write.csv(gee_appendix_table3, row.names = FALSE, quote = FALSE)), collapse = "\n"))

# Convert gee_appendix_table4 to CSV format without index and print
cat(paste0(capture.output(write.csv(gee_appendix_table4, row.names = FALSE, quote = FALSE)), collapse = "\n"))
