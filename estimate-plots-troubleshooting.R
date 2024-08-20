paper_estimates <- read.csv("paper_mx_estimates.csv")

# Plot 1 - MZ & DZ Estimate
library(tidyverse)

e1_mx_data <- paper_estimates %>%
  filter(g == 0) %>%
  dplyr::select("b", "e1", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Estimate = e1)

e2_mx_data <- paper_estimates %>%
  filter(b == 0) %>%
  dplyr::select("g", "e2", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Estimate = e2)

mx_mzdz_data <- rbind(e1_mx_data, e2_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

ggplot(data = mx_mzdz_data, mapping = aes(x = Confounder, y = Estimate, color = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))

# Plot 2 - MZ & DZ vs. DZ-only Estimate
e5_mx_data <- paper_estimates %>%
  filter(g == 0) %>%
  dplyr::select("b", "e5", "PGS") %>%
  mutate(Variable = "CT", Sample = "DZ") %>%
  rename(Confounder = b, Estimate = e5)

e6_mx_data <- paper_estimates %>%
  filter(b == 0) %>%
  dplyr::select("g", "e6", 'PGS') %>%
  mutate(Variable = "SI", Sample = "DZ") %>%
  rename(Confounder = g, Estimate = e6)

mx_dz_data <- rbind(e5_mx_data, e6_mx_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%")))

full_mx_data <- rbind(mx_mzdz_data, mx_dz_data)

ggplot(data = full_mx_data, mapping = aes(x = Confounder, y = Estimate, color = Variable, linetype = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("DZ" = "dotted", "MZ & DZ" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))

# Plot 5
gee_est <- read.csv("paper_gee_estimates.csv")

e1_gee_data <- gee_est %>%
  filter(g == 0) %>%
  dplyr::select("b", "e1", "PGS") %>%
  mutate(Variable = "CT", Sample = "MZ & DZ") %>%
  rename(Confounder = b, Estimate = e1)

e2_gee_data <- gee_est %>%
  filter(b == 0) %>%
  dplyr::select("g", "e2", 'PGS') %>%
  mutate(Variable = "SI", Sample = "MZ & DZ") %>%
  rename(Confounder = g, Estimate = e2)

gee_mzdz_data <- rbind(e1_gee_data, e2_gee_data) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  mutate(Method = "Gee")


mx_mzdz_data <- mx_mzdz_data %>%
  mutate(Method = "OpenMx")

full_method_data <- rbind(mx_mzdz_data, gee_mzdz_data)

ggplot(data = full_method_data, mapping = aes(x = Confounder, y = Estimate, color = Variable, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  facet_wrap(~PGS_percent) +
  scale_linetype_manual(values = c("Gee" = "dotted", "OpenMx" = "solid")) +
  jtools::theme_apa() +
  theme(text = element_text(family = "serif"))
