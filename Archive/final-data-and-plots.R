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

plot1new <- ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif")) +
  scale_color_viridis_d(name = "variable")

plot1 <- ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~PGS_percent) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal(base_family = "CMU Serif")

cairo_pdf("figure3-latex.pdf", width = 6, height = 4)
print(plot1)
dev.off()


ggsave("figure3-latex", width = 6, height = 4, units = "in")

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
plot3 <- ggplot(data_noCT_long, aes(x = SI, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Sibling Interaction",
       y = "Power",
       color = "Parameter", shape = "Parameter") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  theme_minimal(base_family = "CMU Serif")

cairo_pdf("appendix1-latex.pdf", width = 6, height = 4)
print(plot3)
dev.off()

# Plot 4
# Filter the data for b = 0
data_noSI <- ext_power[ext_power$SI== 0, ]

# Gather the data for plotting
data_noSI_long <- gather(data_noSI, key = "variable", value = "value", `CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_noSI_long$variable <- factor(data_noSI_long$variable,
                                  levels = c("CT (m1)", "SI (m2)", "CT (m3)", "SI (m3)"))

# Create the plot
plot4 <- ggplot(data_noSI_long, aes(x = CT, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Cultural Transmission",
       y = "Power",
       color = "Parameter", shape = "Parameter") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  theme_minimal(base_family = "CMU Serif")

cairo_pdf("appendix2-latex.pdf", width = 6, height = 4)
print(plot4)
dev.off()

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
scale_color_manual(values = c("red", "blue")) +
theme_minimal(base_family = "CMU Serif")

cairo_pdf("appendix8-latex.pdf", width = 6, height = 4)
print(plot5)
dev.off()


# Set width and height in inches (1 inch = 25.4 mm)
width <- 200 / 25.4  # Convert 200 mm to inches
height <- 150 / 25.4  # Convert 150 mm to inches

# Export each plot as a PDF with the desired size
ggsave("plot1.pdf", plot = plot1, device = "pdf", width = width, height = height)
ggsave("plot2.pdf", plot = plot2, device = "pdf", width = width, height = height)
ggsave("plot3.pdf", plot = plot3, device = "pdf", width = width, height = height)
ggsave("plot4.pdf", plot = plot4, device = "pdf", width = width, height = height)
ggsave("plot5.pdf", plot = plot5, device = "pdf", width = width, height = height)

ggsave("plot1.jpg", plot = plot1, device = "jpeg", width = width, height = height, dpi = 300)
ggsave("plot2.jpg", plot = plot2, device = "jpeg", width = width, height = height, dpi = 300)
ggsave("plot3.jpg", plot = plot3, device = "jpeg", width = width, height = height, dpi = 300)
ggsave("plot4.jpg", plot = plot4, device = "jpeg", width = width, height = height, dpi = 300)
ggsave("plot5.jpg", plot = plot5, device = "jpeg", width = width, height = height, dpi = 300)


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

# Plots powerpoint
# Plot 3
ext_estimates <- read.csv("2024-06-04_mx_estimates_ext.csv")
ext_power <- read.csv("2024-06-04_mx_power_ext.csv")

# Rename columns with full descriptions
ext_power <- ext_power %>%
  rename(`CT` = g,
         `SI` = b,
         `Cultural Transmission Only` = p1,
         `Sibling Interaction Only` = p2,
         `Cultural Transmission - Combined Model` = p3,
         `Sibling Interaction - Combined Model` = p4)

# Filter the data where CT = 0 for plotting
data_noCT <- ext_power[ext_power$CT == 0, ]

# Gather the data for plotting
data_noCT_long <- gather(data_noCT, key = "variable", value = "value",
                         `Cultural Transmission Only`,
                         `Sibling Interaction Only`,
                         `Cultural Transmission - Combined Model`,
                         `Sibling Interaction - Combined Model`)

# Ensure correct ordering of factor levels
data_noCT_long$variable <- factor(data_noCT_long$variable,
                                  levels = c("Cultural Transmission Only",
                                             "Sibling Interaction Only",
                                             "Cultural Transmission - Combined Model",
                                             "Sibling Interaction - Combined Model"))

# Create the plot for data_noCT
ggplot(data_noCT_long, aes(x = SI, y = value, color = variable, shape = variable)) +
  geom_line(linewidth = 1.2) +  # Adjust line thickness
  geom_point(size = 3) +   # Adjust point size
  labs(x = "Sibling Interaction", y = "Power", color = "Parameter", shape = "Parameter") +  # Ensure both have the same label
  scale_color_manual(values = c("red", "blue", "turquoise", "purple")) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme_minimal()


# Plot 4: Filter the data where SI = 0 for plotting
data_noSI <- ext_power[ext_power$SI == 0, ]

# Gather the data for plotting
data_noSI_long <- gather(data_noSI, key = "variable", value = "value",
                         `Cultural Transmission Only`,
                         `Sibling Interaction Only`,
                         `Cultural Transmission - Combined Model`,
                         `Sibling Interaction - Combined Model`)

# Ensure correct ordering of factor levels
data_noSI_long$variable <- factor(data_noSI_long$variable,
                                  levels = c("Cultural Transmission Only",
                                             "Sibling Interaction Only",
                                             "Cultural Transmission - Combined Model",
                                             "Sibling Interaction - Combined Model"))

# Create the plot
ggplot(data_noSI_long, aes(x = CT, y = value, color = variable, shape = variable)) +
  geom_line(linewidth = 1.2) +  # Adjust line thickness
  geom_point(size = 3) +   # Adjust point size
  labs(x = "Cultural Transmission", y = "Power", color = "Parameter", shape = "Parameter") +  # Ensure both have the same label
  scale_color_manual(values = c("red", "blue", "turquoise", "purple")) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme_minimal()


part1 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p1)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p1),
                color = ifelse(p1 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue", limits = c(0, 1)) +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  ggtitle("CT-Only") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none")

part2 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p2)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p2),
                color = ifelse(p2 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  ggtitle("SI-Only") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none")

part3 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p3)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p3),
                color = ifelse(p3 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  ggtitle("CT - Combined Model") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none")

part4 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p4)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p4),
                color = ifelse(p4 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  ggtitle("SI - Combined Model") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none")

#install.packages("cowplot")
library(cowplot)

# Combine the plots into a grid
combined_plot <- plot_grid(part1, part2, part3, part4, ncol = 2)

# Show the combined plot
print(combined_plot)

# Save the combined plot as a PDF


cairo_pdf("figure2-latex.pdf", width = 6, height = 4)
print(combined_plot)
dev.off()

# Presentation Appendix
plot2 <- ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal(base_family = "CMU Serif")

cairo_pdf("figure4-latex.pdf", width = 6, height = 4)
print(plot2)
dev.off()


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
ggplot(data = full_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal()


filter(full_mzdz_data, Sample == "MZ & DZ", Method == "OpenMx") %>%
  ggplot(data = ., mapping = aes(x = Confounder, y = Power, color = Variable)) +
    geom_line(linewidth = 0.8) +
    geom_point() +
    facet_wrap(~PGS_percent) +
    scale_color_manual(values = c("red", "blue")) +
    theme_minimal()

ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal()


# Vary A

library(OpenMx)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)

# Run function for detailed lots in running text
varya_data <- gnomesims::gnome_mx_simulation(a = sqrt(seq(.4,.8,.1)), ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                             nloci = 100, npgsloci = 10)

# Create seperate data sets for processing
varya_power <- varya_data$power
varya_estimates <- varya_data$params

# Write to .csv files
write.csv(varya_estimates, file = "varya_mx_estimates.csv", row.names = TRUE)
write.csv(varya_power, file = "varya_mx_power.csv", row.names = TRUE)

power_varya <- read.csv("varya_mx_power.csv")

varya1_plot <- power_varya %>%
  ggplot(aes(x = g, y = b, fill = p1)) +
  geom_tile() +
  facet_wrap(~a) +
  geom_text(aes(label = sprintf("%.2f", p1),
                color = ifelse(p1 < 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  theme_minimal(base_family = "CMU Serif")

varya2_plot <-power_varya %>%
  ggplot(aes(x = g, y = b, fill = p2)) +
  geom_tile() +
  facet_wrap(~a) +
  geom_text(aes(label = sprintf("%.2f", p2),
                color = ifelse(p2 < 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  theme_minimal(base_family = "CMU Serif")

varya3_plot <-power_varya %>%
  ggplot(aes(x = g, y = b, fill = p3)) +
  geom_tile() +
  facet_wrap(~a) +
  geom_text(aes(label = sprintf("%.2f", p3),
                color = ifelse(p3 < 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  theme_minimal(base_family = "CMU Serif")

varya4_plot <- power_varya %>%
  ggplot(aes(x = g, y = b, fill = p4)) +
  geom_tile() +
  facet_wrap(~a) +
  geom_text(aes(label = sprintf("%.2f", p4),
                color = ifelse(p4 < 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power") +
  theme_minimal(base_family = "CMU Serif")

cairo_pdf("appendix4-latex.pdf", width = 6, height = 4)
print(varya1_plot)
dev.off()

cairo_pdf("appendix5-latex.pdf", width = 6, height = 4)
print(varya2_plot)
dev.off()

cairo_pdf("appendix6-latex.pdf", width = 6, height = 4)
print(varya3_plot)
dev.off()

cairo_pdf("appendix7-latex.pdf", width = 6, height = 4)
print(varya4_plot)
dev.off()


gnome_effect(g = power_varya$g, b = power_varya$b,
             a = sqrt(power_varya$a), c =  sqrt(power_varya$c),
             e = sqrt(power_varya$e))

# Plotting components
comp_plot_data <- data.frame(setting = 1:4, a = sqrt(seq(.4,1,.2)), c = sqrt(.3), e = sqrt(.3))

# Creating and standardizing the data
comp_plot_data$a_stand <- comp_plot_data$a / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)
comp_plot_data$c_stand <- comp_plot_data$c / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)
comp_plot_data$e_stand <- comp_plot_data$e / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)

# Preparing data for ggplot (pivot for long format)
comp_plot_data_long <- comp_plot_data %>%
  pivot_longer(cols = c(a, c, a_stand, c_stand),
               names_to = c(".value", "type"),
               names_pattern = "(.*)(_stand)?")

# Plot unstandardized
p1 <- ggplot(comp_plot_data_long, aes(x = setting)) +
  geom_line(aes(y = a, color = "blue")) +
  geom_line(aes(y = c, color = "red")) +
  # Add white circles behind the text
  geom_point(aes(y = a), color = "white", size = 6) +  # White background for letter "a"
  geom_point(aes(y = c), color = "white", size = 6) +  # White background for letter "c"
  geom_text(aes(y = a, label = "a", color = "blue"), size = 4, family = "CMU Serif") +   # Use "a" as the point
  geom_text(aes(y = c, label = "c", color = "red"), size = 4, family = "CMU Serif") +    # Use "c" as the point
  labs(x = "Setting", y = "Unstandardized Values") +
  scale_color_identity(guide = "none") +
  ylim(0,1) +
  theme_minimal(base_family = "CMU Serif")

# Plot standardized
p2 <- ggplot(comp_plot_data_long, aes(x = setting)) +
  geom_line(aes(y = a_stand, color = "blue")) +
  geom_line(aes(y = c_stand, color = "red")) +
  # Add white circles behind the text
  geom_point(aes(y = a_stand), color = "white", size = 6) +  # White background for letter "a"
  geom_point(aes(y = c_stand), color = "white", size = 6) +  # White background for letter "c"
  geom_text(aes(y = a_stand, label = "a", color = "blue", family = "CMU Serif"), size = 4) +   # Use "a" as the point
  geom_text(aes(y = c_stand, label = "c", color = "red", family = "CMU Serif"), size = 4) +    # Use "c" as the point
  labs(x = "Setting", y = "Standardized Values") +
  scale_color_identity(guide = "none") +
  ylim(0,1) +
  theme_minimal(base_family = "CMU Serif")

# Display the plots side by side
library(gridExtra)

cairo_pdf("appendix3-latex.pdf", width = 6, height = 2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# Plots for parameter settings
par(mfrow = c(1, 2))
plot(comp_plot_data$setting, comp_plot_data$a, type = "b", pch = "a",
     col = "blue", ylim = c(0, 1), xlab = "Setting",
     ylab = "Unstandardized Values", xaxt = "n")
axis(1, at = comp_plot_data$setting)
lines(comp_plot_data$setting, comp_plot_data$c, type = "b", pch = "c", col = "red")
plot(comp_plot_data$setting, comp_plot_data$a_stand, type = "b", pch = "a",
     col = "blue", ylim = c(0, 1), xlab = "Setting",
     ylab = "Standardized Values", xaxt = "n")
lines(comp_plot_data$setting, comp_plot_data$c_stand, type = "b", pch = "c", col = "red")
axis(1, at = comp_plot_data$setting)
par(mfrow = c(1, 1))

# Same in ggplot
# Reshape data for plotting
long_data <- comp_plot_data %>%
  pivot_longer(cols = c(a, c, a_stand, c_stand),
               names_to = c(".value", "type"),
               names_pattern = "([a-zA-Z]+)_(.*)")

# Create the ggplot with facets
ggplot(long_data, aes(x = setting, y = value, color = type, group = type)) +
  geom_line(aes(linetype = type), size = 1) +  # Line for a and c
  geom_point(aes(shape = type), size = 3) +   # Points for a and c
  facet_wrap(~type, scales = "free_y", labeller = labeller(type = c("a" = "Unstandardized Values", "c" = "Standardized Values"))) +
  labs(x = "Setting", y = "Values") +
  scale_color_manual(values = c("blue", "red")) +  # Colors for a and c
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if desired

# Appendix Tables
ext_estimates <- read.csv("2024-06-04_mx_estimates_ext.csv")
ext_power <- read.csv("2024-06-04_mx_power_ext.csv")

table1_appendix <- ext_power %>%
  select(g, b, p1:p4, p13:p16, Smz, Sdz) %>%
  mutate(Smz = scales::percent(Smz,accuracy=1),
         Sdz = scales::percent(Sdz,accuracy=1))

cat(paste0(capture.output(write.csv(table1_appendix, row.names = FALSE, quote = FALSE)), collapse = "\n"))

table2_appendix <- ext_estimates %>%
  select(g, b, e1:e4, e13:e16, Smz, Sdz) %>%
  mutate(Smz = scales::percent(Smz,accuracy=1),
         Sdz = scales::percent(Sdz,accuracy=1))

cat(paste0(capture.output(write.csv(table2_appendix, row.names = FALSE, quote = FALSE)), collapse = "\n"))
