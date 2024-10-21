# Install all packages
# detach("package:gnomesims", unload = TRUE)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)
library(OpenMx)

# Run function for detailed lots in running text
# paper_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02), nloci = 100, npgsloci = c(2, 5, 10, 15))
# gee_data <- gnomesims::gnome_gee_simulation(ct = seq(0,.1,.02), si = seq(0,.1,.02), nloci = 100, npgsloci = c(2, 5, 10, 15), cmethod = "exchangeable")

# Create seperate data sets for processing
# paper_power <- paper_data$power
# paper_estimates <- paper_data$params
# gee_power <- gee_data$power
# gee_estimates <- gee_data$params

# Write to .csv files
# write.csv(paper_estimates, file = "paper_mx_estimates.csv", row.names = TRUE)
# write.csv(paper_power, file = "paper_mx_power.csv", row.names = TRUE)
# write.csv(gee_estimates, file = "paper_gee_estimates.csv", row.names = TRUE)
# write.csv(gee_power, file = "paper_gee_power.csv", row.names = TRUE)

# Start here if files are already created
paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")

gee_power <- read.csv("paper_gee_power.csv")

# Run function for detailed lots in running text
# varya_data <- gnomesims::gnome_mx_simulation(a = sqrt(seq(.4,.8,.1)), ct = seq(0,.1,.02), si = seq(0,.1,.02), nloci = 100, npgsloci = 10)

# Create seperate data sets for processing
# varya_power <- varya_data$power
# varya_estimates <- varya_data$params

# Write to .csv files
# write.csv(varya_estimates, file = "varya_mx_estimates.csv", row.names = TRUE)
# write.csv(varya_power, file = "varya_mx_power.csv", row.names = TRUE)

# Start here if file is already created
power_varya <- read.csv("varya_mx_power.csv")
