# Install all packages

devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)
library(OpenMx)

# Run function for detailed lots in running text
paper_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.02),
                  si = seq(0,.1,.02), nloci = 100, npgsloci = c(2, 5, 10, 15))

gee_data <- gnomesims::gnome_gee_simulation(ct = seq(0,.1,.02),
                  si = seq(0,.1,.02), nloci = 100, npgsloci = c(2, 5, 10, 15),
                  cmethod = "exchangeable")

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

# Revision - Create data for different sample sizes
# Define sample sizes
sample_sizes <- c(seq(10000, 50000, 10000))
sample_sizes_dz <- sample_sizes/2
unique_sample_sizes <- unique(c(sample_sizes, sample_sizes_dz))

sample_list <- list()
for (n in unique_sample_sizes) {
   if (length(sample_list) < which(unique_sample_sizes == n)) {
     result <- gnome_mx_simulation(nmz = n, ndz = n, ct = 0.03, si = 0.03, nloci = 100, npgsloci = 10)
     sample_list <- append(sample_list, list(result))
   }
 }

# Extract power and parameters
sample_power <- do.call(rbind, lapply(sample_list, `[[`, "power"))
sample_params <- do.call(rbind, lapply(sample_list, `[[`, "params"))

write.csv(sample_power, "sample_power.csv", row.names = TRUE)
write.csv(sample_params, "sample_params.csv", row.names = TRUE)


# Start here if files are already created
paper_power <- read.csv("paper_mx_power.csv")
paper_estimates <- read.csv("paper_mx_estimates.csv")
