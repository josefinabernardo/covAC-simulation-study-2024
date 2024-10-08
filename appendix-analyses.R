library(OpenMx)
devtools::install_github("josefinabernardo/gnomesims", force = TRUE)
library(gnomesims)

# Results Sibling Mean
gnomesims::gnome_mx_simulation(a = sqrt(.74), c = sqrt(.08), e = sqrt(.18),
                               nmz = 2451, ndz = 4569, nloci = 100,
                               npgsloci = 12, ct = seq(.1,.13,.01), si = seq(.1,.13,.01))

