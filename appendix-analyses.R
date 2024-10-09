library(OpenMx)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)

# Results Sibling Mean
dezeeuw_data <- gnomesims::gnome_mx_simulation(a = sqrt(.74), c = sqrt(.08),
                      e = sqrt(.18), nmz = 2451, ndz = 4569, nloci = 100,
                      npgsloci = 12, ct = seq(0,.12,.01), si = seq(0,.12,.01))

# Create seperate data sets for processing
dezeeuw_power <- dezeeuw_data$power
dezeeuw_estimates <- dezeeuw_data$params

# Write to .csv files
write.csv(dezeeuw_estimates, file = "dezeeuw_estimates.csv", row.names = TRUE)
write.csv(dezeeuw_power, file = "dezeeuw_power.csv", row.names = TRUE)

dezeeuw_power <- read.csv("dezeeuw_power.csv")

dezeeuw_power %>%
  filter(g == b) %>%
  select(g, b, p1, p2, p3, p4, Smz, Sdz)


