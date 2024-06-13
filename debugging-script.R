# DEBUGGING
data_list <- dolan_simulation_function(a = sqrt(.4), c = sqrt(.3),
                                       e = sqrt(.3), nloci = 100, npgsloci = c(2,5,10,15),
                                       ct = sqrt(c(0, .0025, 0.01)), si = sqrt(c(0, .0025, 0.01)))


# Extract data sets
mx_estimates <- drop_na(data_list[[1]])
mx_power <- drop_na(data_list[[2]])
gee_estimates <- drop_na(data_list[[3]])
gee_power <- drop_na(data_list[[4]])

# Re-name columns
setnames = c('nmz','ndz','a','c','e','g','b','x','PGS','A')
colnames(mx_estimates) <- c(setnames, paste0("e", 1:16))
colnames(mx_power) <- c(setnames, paste0("p", 1:16))
colnames(gee_estimates) <- c(setnames, paste0("e", 1:9))
colnames(gee_power) <- c(setnames, paste0("p", 1:9))

# Use effect size function on the data sets
mxestimates <- mx_estimates %>%
  mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
         Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)
mxpower <- mx_power %>%
  mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
         Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)
geeestimates <- gee_estimates %>%
  mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
         Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)
geepower <- gee_power %>%
  mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
         Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)

# Write data frames to CSV files
write.csv(mxestimates, file = "debug_mx_estimates.csv", row.names = TRUE)
write.csv(mxpower, file = "debug_mx_power.csv", row.names = TRUE)
write.csv(geeestimates, file = "debug_gee_estimates.csv", row.names = TRUE)
write.csv(geepower, file = "debug_gee_power.csv", row.names = TRUE)

# Debugged data
app_mx_estimates <- read.csv("debug_mx_estimates.csv")
app_mx_power <- read.csv("debug_mx_power.csv")
app_gee_estimates <- read.csv("debug_gee_estimates.csv")
app_gee_power <- read.csv("debug_gee_power.csv")

library(tidyverse)

# PLOTTING

# Start with some overview plots

# MZ DZ (Model 1 vs Model 3)
app_mx_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1, `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4) %>%
  filter(CT == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = SI, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "SI",
       y = "Value") +
  theme_minimal()

app_mx_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1, `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4) %>%
  filter(SI == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "CT",
       y = "Value") +
  theme_minimal()

# SOMETHING IS WRONG HERE... is it the /2? CT model 1 is too low
app_gee_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p6, `SI (m2)` = p7, `CT (m3)` = p8, `SI (m3)` = p9) %>%
  filter(CT == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = SI, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "SI",
       y = "Value") +
  theme_minimal()

app_gee_power %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p6, `SI (m2)` = p7, `CT (m3)` = p8, `SI (m3)` = p9) %>%
  filter(SI == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "CT",
       y = "Value") +
  theme_minimal()

# Same for estimates
app_mx_estimates %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e1, `SI (m2)` = e2, `CT (m3)` = e3, `SI (m3)` = e4) %>%
  filter(CT == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = SI, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "SI",
       y = "Value") +
  theme_minimal()

app_mx_estimates %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e1, `SI (m2)` = e2, `CT (m3)` = e3, `SI (m3)` = e4) %>%
  filter(SI == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "CT",
       y = "Value") +
  theme_minimal()

app_gee_estimates %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e6, `SI (m2)` = e7, `CT (m3)` = e8, `SI (m3)` = e9) %>%
  filter(CT == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = SI, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "SI",
       y = "Value") +
  theme_minimal()

app_gee_estimates %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = e6, `SI (m2)` = e7, `CT (m3)` = e8, `SI (m3)` = e9) %>%
  filter(SI == 0) %>%
  mutate(PGS_percent = factor(scales::percent(PGS), levels = c("2%", "5%", "10%", "15%"))) %>%
  gather(key = "variable", value = "value",`CT (m1)`, `SI (m2)`, `CT (m3)`, `SI (m3)`) %>%
  ggplot(aes(x = CT, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~PGS_percent) +
  labs(title = NULL,
       x = "CT",
       y = "Value") +
  theme_minimal()

# Comparing combined sample (Model 1 OpenMx vs. Model 3 Gee)
mzdz_power_mx <- app_mx_power %>%
  select(b, g, p1, p2, p3 , p4) %>%
  rename(`CT` = g, `SI` = b, `CT (m1)` = p1,
         `SI (m2)` = p2, `CT (m3)` = p3, `SI (m3)` = p4)

mzdz_power_gee <- app_gee_power

mzdz_est_mx <- app_mx_estimates

mzdz_est_gee <- app_gee_estimates


# Comparing DZ-only sample (Model 4 OpenMx vs. Model 2 Gee)
dz_power_mx <- app_mx_power

dz_power_gee <- app_gee_power

dz_est_mx <- app_mx_estimates

dz_est_gee <- app_gee_estimates
