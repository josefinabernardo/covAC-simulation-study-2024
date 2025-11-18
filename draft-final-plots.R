# Install all packages
devtools::install_github("josefinabernardo/gnomesims")
library(gnomesims)
library(OpenMx)

# Run function for detailed lots in running text
paper_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                        nloci = 100,
                                        npgsloci = c(2, 5, 10, 15))

# Create seperate data sets for processing
paper_power <- paper_data$power
paper_estimates <- paper_data$params

# Write to .csv files
write.csv(paper_estimates, file = "paper_mx_estimates.csv", row.names = TRUE)
write.csv(paper_power, file = "paper_mx_power.csv", row.names = TRUE)

# Load in the data sets
paper_power <- read.csv("2024_08_14_mx_power.csv")
paper_estimates <- read.csv("2024_08_14_mx_estimates.csv")

#paper_power <- read.csv("paper_mx_power.csv")
#paper_estimates <- read.csv("paper_mx_estimates.csv")

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
  dplyr::select("b", "p13", "PGS") %>%
  mutate(Variable = "CT", Sample = "DZ") %>%
  rename(Confounder = b, Power = p13)

p6_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p14", 'PGS') %>%
  mutate(Variable = "SI", Sample = "DZ") %>%
  rename(Confounder = g, Power = p14)

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
