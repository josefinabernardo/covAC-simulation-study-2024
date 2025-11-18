library(tidyverse)
library(scales)
library(gnomesims)

# Appendix data
app_mx_estimates <- read.csv("2024-05-08_mx_estimates_appendix.csv")
app_mx_power <- read.csv("2024-05-08_mx_power_appendix.csv")
app_gee_estimates <- read.csv("2024-05-08_gee_estimates_appendix.csv")
app_gee_power <- read.csv("2024-05-08_gee_power_appendix.csv")


app_mx_estimates <- app_mx_estimates %>%
  select(-X) %>%
  mutate(Smz = percent(round(gnome_effect(a, c, e, g, b)$mz, 4)),
         Sdz = percent(round(gnome_effect(a, c, e, g, b)$dz, 4)))

app_mx_power <- app_mx_power %>%
  select(-X) %>%
  mutate(Smz = percent(round(gnome_effect(a, c, e, g, b)$mz, 4)),
         Sdz = percent(round(gnome_effect(a, c, e, g, b)$dz, 4)))

app_gee_estimates <- app_gee_estimates %>%
  select(-X) %>%
  mutate(Smz = percent(round(gnome_effect(a, c, e, g, b)$mz, 4)),
         Sdz = percent(round(gnome_effect(a, c, e, g, b)$dz, 4)))

app_gee_power <- app_gee_power %>%
  select(-X) %>%
  mutate(Smz = percent(round(gnome_effect(a, c, e, g, b)$mz, 4)),
         Sdz = percent(round(gnome_effect(a, c, e, g, b)$dz, 4)))

# Write to files
write.csv(app_mx_estimates, file = "2024-05-13_app_mx_estimates.csv")
write.csv(app_mx_power, file = "2024-05-13_app_mx_power.csv")
write.csv(app_gee_estimates, file = "2024-05-13_app_gee_estimates.csv")
write.csv(app_gee_power, file = "2024-05-13_app_gee_power.csv")

