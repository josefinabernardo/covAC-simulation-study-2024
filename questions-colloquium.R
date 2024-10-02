# Load package
library(OpenMx)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)

# Run function for detailed lots in running text
am_data <- gnomesims::gnome_mx_simulation(ct = seq(0,.1,.05), si = seq(0,.1,.05),
                                             nloci = 100,
                                             npgsloci = 10,
                                             assortm = seq(0,.4,.1))
