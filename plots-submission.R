# Packages
library(tidyverse)
library(cowplot)
library(gridExtra)

# Load in data sets
paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")

# RUNNING TEXT
# Plot 2 - Statistical Power Relative to effect sizes for cultural transmission and sibling interaction
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

# Combine the plots into a grid
plot2 <- plot_grid(part1, part2, part3, part4, ncol = 2)

# Plot 3 - Type I Error Rate for Four Different Strengths of PGS Predictive Power in the MZ & DZ Sample
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

plot3 <- ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~PGS_percent) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal(base_family = "CMU Serif")

# Plot 4 - Type I Error Rate for Four Different Strengths of PGS Predictive Power in the MZ & DZ vs. DZ-only Sample
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

plot4 <- ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal(base_family = "CMU Serif")

# APPENDIX
# Plot 1
appendix1 <- ggplot(data_noCT_long, aes(x = SI, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Sibling Interaction",
       y = "Power",
       color = "Parameter", shape = "Parameter") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  theme_minimal(base_family = "CMU Serif")

# Plot 2
appendix2 <- ggplot(data_noSI_long, aes(x = CT, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Cultural Transmission",
       y = "Power",
       color = "Parameter", shape = "Parameter") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  theme_minimal(base_family = "CMU Serif")

# Plot 3
# Creating and standardizing the data
comp_plot_data <- data.frame(setting = 1:4, a = sqrt(seq(.4,1,.2)), c = sqrt(.3), e = sqrt(.3))
comp_plot_data$a_stand <- comp_plot_data$a / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)
comp_plot_data$c_stand <- comp_plot_data$c / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)
comp_plot_data$e_stand <- comp_plot_data$e / (comp_plot_data$a + comp_plot_data$c + comp_plot_data$e)

# Preparing data for ggplot (pivot for long format)
comp_plot_data_long <- comp_plot_data %>%
  pivot_longer(cols = c(a, c, a_stand, c_stand),
               names_to = c(".value", "type"),
               names_pattern = "(.*)(_stand)?")

# Plot unstandardized
app_part1 <- ggplot(comp_plot_data_long, aes(x = setting)) +
  geom_line(aes(y = a, color = "blue")) +
  geom_line(aes(y = c, color = "red")) +
  geom_point(aes(y = a), color = "white", size = 6) +  # White background for letter "a"
  geom_point(aes(y = c), color = "white", size = 6) +  # White background for letter "c"
  geom_text(aes(y = a, label = "a", color = "blue"), size = 4, family = "CMU Serif") +   # Use "a" as the point
  geom_text(aes(y = c, label = "c", color = "red"), size = 4, family = "CMU Serif") +    # Use "c" as the point
  labs(x = "Setting", y = "Unstandardized Values") +
  scale_color_identity(guide = "none") +
  ylim(0,1) +
  theme_minimal(base_family = "CMU Serif")

# Plot standardized
app_part2 <- ggplot(comp_plot_data_long, aes(x = setting)) +
  geom_line(aes(y = a_stand, color = "blue")) +
  geom_line(aes(y = c_stand, color = "red")) +
  geom_point(aes(y = a_stand), color = "white", size = 6) +  # White background for letter "a"
  geom_point(aes(y = c_stand), color = "white", size = 6) +  # White background for letter "c"
  geom_text(aes(y = a_stand, label = "a", color = "blue", family = "CMU Serif"), size = 4) +   # Use "a" as the point
  geom_text(aes(y = c_stand, label = "c", color = "red", family = "CMU Serif"), size = 4) +    # Use "c" as the point
  labs(x = "Setting", y = "Standardized Values") +
  scale_color_identity(guide = "none") +
  ylim(0,1) +
  theme_minimal(base_family = "CMU Serif")

grid.arrange(app_part1, app_part2, ncol = 2)

# Plot 4
appendix4 <- power_varya %>%
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

# Plot 5
appendix5 <-power_varya %>%
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

# Plot 6
appendix6 <-power_varya %>%
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

# Plot 7
appendix7 <- power_varya %>%
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

# Plot 8
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

appendix8 <- ggplot(data = full_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal(base_family = "CMU Serif")
