# Packages
library(tidyverse)
library(cowplot)
library(gridExtra)
library(extrafont)
library(Cairo)

# Load in data sets
paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")

sample_power <- read.csv("sample_power.csv")
sample_params <- read.csv("sample_params.csv")


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
  labs(x = NULL, y = NULL, fill = "Power") +
  ggtitle("Power Cultural Transmission\n(Using Parent PGS-Only)") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(margin = margin(t = -5)),
        axis.text.y = element_text(margin = margin(r = -5)))

part2 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p2)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p2),
                color = ifelse(p2 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = NULL, y = NULL, fill = "Power") +
  ggtitle("Power Sibling Interaction\n(Using Sibling PGS-Only)") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(margin = margin(t = -5)),
        axis.text.y = element_text(margin = margin(r = -5)))

part3 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p3)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p3),
                color = ifelse(p3 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = NULL, y = NULL, fill = "Power") +
  ggtitle("Power Cultural Transmission\n(Using Sibling & Parent PGS)") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(margin = margin(t = -5)),
        axis.text.y = element_text(margin = margin(r = -5)))

part4 <- paper_power %>%
  filter(PGS == 0.1) %>%
  ggplot(aes(x = g, y = b, fill = p4)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p4),
                color = ifelse(p4 > 0.5, "white", "black")),  # Dynamic text color
            size = 3, family = "CMU Serif") +
  scale_fill_gradient(low = "red", high = "blue") +  # Two-colored gradient
  scale_color_identity() +  # Use identity scale for text color (white/black)
  labs(x = NULL, y = NULL, fill = "Power") +
  ggtitle("Power Sibling Interaction\n(Using Sibling & Parent PGS)") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(margin = margin(t = -5)),
        axis.text.y = element_text(margin = margin(r = -5)))

# Combine the plots into a grid
plot2_int_1 <- plot_grid(part1, part2, ncol = 2, align = "hv")
plot2_int_2 <- plot_grid(part3, part4, ncol = 2, align = "hv")
plot2_int <- plot_grid(part1, part2, part3, part4, ncol = 2, align = "hv")

Figure2 <- ggdraw() +
  draw_plot(plot2_int_1, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_label("Cultural Transmission", x = 0.5, y = 0, vjust = -0.5, fontfamily = "CMU Serif") +
  draw_label("Sibling Interaction", x = 0, y = 0.5, angle = 90, vjust = 1.5, fontfamily = "CMU Serif")
Figure2

Figure3 <- ggdraw() +
  draw_plot(plot2_int_2, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_label("Cultural Transmission", x = 0.5, y = 0, vjust = -0.5, fontfamily = "CMU Serif") +
  draw_label("Sibling Interaction", x = 0, y = 0.5, angle = 90, vjust = 1.5, fontfamily = "CMU Serif")
Figure3

#  labs(x = "Cultural Transmission", y = "Sibling Interaction", fill = "Power")

# Plot 3 - Type I Error Rate for Four Different Strengths of PGS Predictive Power in the MZ & DZ Sample
p1_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p1", "PGS") %>%
  mutate(Variable = "Cultural Transmission", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p1)

p2_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p2", 'PGS') %>%
  mutate(Variable = "Sibling Interaction", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Power = p2)

mx_mzdz_data <- rbind(p1_mx_data, p2_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

axis_lines_mzdz <- mx_mzdz_data %>%
  group_by(PGS_percent) %>%
  summarise(x_start = 0, x_end = max(Confounder, na.rm = TRUE),
            y_start = 0, y_end = 0.7) %>%
  ungroup()

Figure4 <- ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Power, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~PGS_percent) +
  geom_segment(data = axis_lines_mzdz,
               aes(x = x_start, xend = x_end, y = 0, yend = 0),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 10),
               color = "gray50", linewidth = 0.5) +
  geom_segment(data = axis_lines_mzdz,
               aes(x = 0, xend = 0, y = y_start, yend = y_end),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 10),
               color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = c("red", "blue")) +
  labs(color = "Modeled Source of\nAC Covariance") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
Figure4

# Plot 5 - Type I Error Rate for Four Different Strengths of PGS Predictive Power in the MZ & DZ vs. DZ-only Sample
p3_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p3", "PGS") %>%
  mutate(Variable = "Cultural Transmission", Sample = "MZ & DZ") %>%
  rename(Parameter = g, Power = p3)

p4_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p4", 'PGS') %>%
  mutate(Variable = "Sibling Interaction", Sample = "MZ & DZ") %>%
  rename(Parameter = b, Power = p4)

p7_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p7", "PGS") %>%
  mutate(Variable = "Cultural Transmission", Sample = "DZ") %>%
  rename(Parameter = g, Power = p7)

p8_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p8", 'PGS') %>%
  mutate(Variable = "Sibling Interaction", Sample = "DZ") %>%
  rename(Parameter = b, Power = p8)

full_mx_data <- rbind(p3_mx_data, p4_mx_data, p7_mx_data, p8_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

axis_lines_pgs <- full_mx_data %>%
  group_by(PGS_percent) %>%
  summarise(x_start = 0, x_end = max(Parameter),
    y_start = 0, y_end = 0.8) %>%
  ungroup()

Figure5 <- ggplot(data = full_mx_data) +
  geom_line(linewidth = 0.8, mapping = aes(x = Parameter, y = Power, color = Variable, linetype = Sample)) +
  geom_point(size = 1.5, mapping = aes(x = Parameter, y = Power, color = Variable)) +
  facet_wrap(~PGS_percent) +
  geom_segment(data = axis_lines_pgs,
               aes(x = x_start, xend = x_end, y = 0, yend = 0),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 10),
               color = "gray50", linewidth = 0.5) +
  geom_segment(data = axis_lines_pgs,
               aes(x = 0, xend = 0, y = y_start, yend = y_end),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 10),
               color = "gray50", linewidth = 0.5) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  scale_color_manual(values = c("red", "blue")) +
  labs(color = "Modeled Source of\nAC Covariance") +
  theme_minimal(base_family = "CMU Serif") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
Figure5

# Prepare sample data
sample_data <- sample_power %>%
  select(nmz, ndz, p3, p4, p7, p8) %>%
  mutate(N = nmz + ndz) %>%
  rename(c("Sample Size" = "N", "Cultural Transmission MZ & DZ" = "p3",
           "Sibling Interaction MZ & DZ" = "p4",
           "Cultural Transmission DZ" = "p7",
           "Sibling Interaction DZ" = "p8"))

sample_data_long <- sample_data %>%
  mutate(Sample_Size_Original = `Sample Size`) %>%
  pivot_longer(cols = c("Cultural Transmission MZ & DZ",
                        "Sibling Interaction MZ & DZ",
                        "Cultural Transmission DZ",
                        "Sibling Interaction DZ"),
               names_to = "variable", values_to = "value") %>%
  mutate(`Sample Size` = ifelse(str_detect(variable, "DZ$") & !str_detect(variable, "MZ & DZ"),
                                Sample_Size_Original / 2,
                                Sample_Size_Original)) %>%
  select(-Sample_Size_Original)

common_sample_sizes <- sample_data_long %>%
  group_by(`Sample Size`) %>%
  summarise(has_MZ_DZ = any(str_detect(variable, "MZ & DZ")),
            has_DZ = any(str_detect(variable, "DZ$") & !str_detect(variable, "MZ & DZ"))) %>%
  filter(has_MZ_DZ & has_DZ) %>%
  pull(`Sample Size`)

sample_data_long <- sample_data_long %>%
  filter(`Sample Size` %in% common_sample_sizes) %>%
  filter(`Sample Size` %in% c(10000, 20000, 30000, 40000, 50000))

Figure6 <- ggplot(sample_data_long, aes(x = `Sample Size`, y = value,
                             color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Sample Size",
       y = "Power",
       color = "Modeled Source of\nAC Covariance",
       shape = "Modeled Source of\nAC Covariance") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  scale_x_continuous(limits = c(5000, 55000),
                     breaks = seq(0, 55000, by = 10000),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  theme_minimal(base_family = "CMU Serif") +
  theme(
    axis.line = element_line(color = "gray50", linewidth = 0.5),
    axis.ticks = element_line(color = "gray50"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    axis.line.x =
      element_line(arrow = arrow(length = unit(0.3, "cm"),
                                 angle = 10,
                                 type="closed")),
    axis.line.y =
      element_line(arrow = arrow(length = unit(0.3, "cm"),
                                 angle = 10,
                                 type="closed")),
   legend.text = element_text(size = 11),
   legend.title = element_text(size = 13)
  )
Figure6

# Export
# Save all plots as .jpgs
jpeg("Figure2.jpg", width = 6, height = 4, units = "in", res = 600)
print(plot2)
dev.off()

jpeg("Figure2_1.jpg", width = 6, height = 2.5, units = "in", res = 600)
print(plot2_1)
dev.off()

jpeg("Figure2_2.jpg", width = 6, height = 2.5, units = "in", res = 600)
print(plot2_2)
dev.off()

jpeg("Figure3.jpg", width = 8, height = 5, units = "in", res = 600)
print(plot3)
dev.off()

jpeg("Figure4.jpg", width = 8, height = 5, units = "in", res = 300)
print(plot4)
dev.off()

jpeg("Figure6.jpg", width = 8, height = 4, units = "in", res = 300)
print(plot6)
dev.off()

# APPENDIX
ext_power <- read.csv("2024-06-04_mx_power_ext.csv")

ext_power <- ext_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1,
         `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4)

data_noCT <- ext_power[ext_power$CT == 0, ]

# Gather the data for plotting
data_noCT_long <- gather(data_noCT, key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_noCT_long$variable <- factor(data_noCT_long$variable,
                                  levels = c("CT (m1)", "SI (m2)", "CT (m3)", "SI (m3)"))


# Plot 1
appendix1 <- ggplot(data_noCT_long, aes(x = SI, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Sibling Interaction",
       y = "Power",
       color = "Modeled Source of\nAC Covariance", shape = "Modeled Source of\nAC Covariance") +
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(4)) +
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  theme_minimal(base_family = "CMU Serif")

# Plot 2
data_noSI <- ext_power[ext_power$SI== 0, ]

# Gather the data for plotting
data_noSI_long <- gather(data_noSI, key = "variable", value = "value", `CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`)
data_noSI_long$variable <- factor(data_noSI_long$variable,
                                  levels = c("CT (m1)", "SI (m2)", "CT (m3)", "SI (m3)"))

appendix2 <- ggplot(data_noSI_long, aes(x = CT, y = value, color = variable, shape = variable)) +
  geom_line() +
  geom_point() +
  labs(title = NULL,
       x = "Cultural Transmission",
       y = "Power",
       color = "Modeled Source of\nAC Covariance", shape = "Modeled Source of\nAC Covariance") +
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

appendix3 <- plot_grid(app_part1, app_part2, ncol = 2)

# Plot 4
power_varya <- read.csv("varya_mx_power.csv")

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
appendix5 <- power_varya %>%
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
appendix6 <- power_varya %>%
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
gee_power <- read.csv("paper_gee_power.csv")

p1_gee_data <- gee_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p1", "PGS") %>%
  mutate(Variable = "Cultural Transmission", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Power = p1)

p2_gee_data <- gee_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p2", 'PGS') %>%
  mutate(Variable = "Sibling Interaction", Sample = "MZ & DZ") %>%
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
  labs(color = "Modeled Source of\nAC Covariance") +
  theme_minimal(base_family = "CMU Serif")

# Print all plots
print(plot2) # 4-in-1 Power
print(plot3) # Faceted PGS
print(plot4) # Faceted PGS

print(appendix1) # Line graph
print(appendix2) # Line graph
print(appendix3) # Variance components
print(appendix4) # 4-in-1 Power
print(appendix5)
print(appendix6)
print(appendix7)
print(appendix8) #


# Save all plots as .jpgs
jpeg("Figure2.jpg", width = 6, height = 4, units = "in", res = 600)
print(plot2)
dev.off()

jpeg("Figure2_1.jpg", width = 6, height = 2.5, units = "in", res = 600)
print(plot2_1)
dev.off()

jpeg("Figure2_2.jpg", width = 6, height = 2.5, units = "in", res = 600)
print(plot2_2)
dev.off()

jpeg("Figure3.jpg", width = 8, height = 5, units = "in", res = 600)
print(plot3)
dev.off()

jpeg("Figure4.jpg", width = 8, height = 5, units = "in", res = 300)
print(plot4)
dev.off()

jpeg("Figure6.jpg", width = 8, height = 4, units = "in", res = 300)
print(plot6)
dev.off()

jpeg("appendix1-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix1)
dev.off()

jpeg("appendix2-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix2)
dev.off()

jpeg("appendix3-submission.jpg", width = 6, height = 3, units = "in", res = 300)
print(appendix3)
dev.off()

jpeg("appendix4-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix4)
dev.off()

jpeg("appendix5-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix5)
dev.off()

jpeg("appendix6-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix6)
dev.off()

jpeg("appendix7-submission.jpg", width = 6, height = 4, units = "in", res = 300)
print(appendix7)
dev.off()

jpeg("appendix8-submission.jpg", width = 7, height = 5, units = "in", res = 300)
print(appendix8)
dev.off()

# Save all plots as .eps

# Packages
#install.packages("extrafont")
library(extrafont)
font_import(prompt = FALSE)
loadfonts(device = "postscript")

# Save plot2 as .eps
postscript("Figure2.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(plot2)
dev.off()

# Save plot3 as .eps
postscript("Figure3.eps", width = 7, height = 5, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(plot3)
dev.off()

# Save plot4 as .eps
postscript("Figure4.eps", width = 7, height = 5, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(plot4)
dev.off()

# Save appendix1 as .eps
postscript("Figure5.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix1)
dev.off()

# Save appendix2 as .eps
postscript("Figure6.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix2)
dev.off()

# Save appendix3 as .eps
postscript("Figure7.eps", width = 6, height = 3, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix3)
dev.off()

# Save appendix4 as .eps
postscript("Figure8.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix4)
dev.off()

# Save appendix5 as .eps
postscript("Figure9.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix5)
dev.off()

# Save appendix6 as .eps
postscript("Figure10.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix6)
dev.off()

# Save appendix7 as .eps
postscript("Figure11.eps", width = 6, height = 4, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix7)
dev.off()

# Save appendix8 as .eps
postscript("Figure6.eps", width = 7, height = 5, onefile = FALSE, paper = "special", horizontal = FALSE, family = "CMU Serif")
print(appendix8)
dev.off()



