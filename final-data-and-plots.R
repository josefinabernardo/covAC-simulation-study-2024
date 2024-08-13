# Run function for detailed lots in running text
paper_data <- dolan_simulation_function(ct = seq(0,.1,.02), si = seq(0,.1,.02),
                                  nloci = 100,
                                  npgsloci = c(2, 5, 10, 15))

# Create seperate data sets for processing
paper_power <- paper_data$power
paper_estimates <- paper_data$params

paper_power_filtered <- paper_power %>%
  filter(g %in% seq(0,.1,.02) & b %in% seq(0,.1,.02))

paper_estimates_filtered <- paper_estimates %>%
  filter(g %in% seq(0,.1,.02)& b %in% seq(0,.1,.02))

# Write to .csv files
write.csv(paper_estimates, file = "paper_mx_estimates.csv", row.names = TRUE)
write.csv(paper_power, file = "paper_mx_power.csv", row.names = TRUE)

# Load in the data sets
paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")

# Plot 1 - MZ & DZ Power
library(tidyverse)

p1_mx_data <- paper_power %>%
  filter(g == 0) %>%
  dplyr::select("b", "p1", "PGS") %>%
  mutate(Variable = "CT")

p2_mx_data <- paper_power %>%
  filter(b == 0) %>%
  dplyr::select("g", "p2", 'PGS') %>%
  mutate(Variable = "SI")
