# First R coding from scratch

2 + 3 

anna <- 2 + 3 # assign an operation to an object
chiara <- 4 + 6 

anna + chiara

filippo <- c(0.2, 0.4, 0.6, 0.8, 0.9)     # array
luca <- c(100, 80, 60, 50, 10)

plot(luca, filippo)
plot(luca, filippo, pch=19)
plot(luca, filippo, pch=19, col="blue")
plot(luca, filippo, pch=19, col="blue", cex=2)
plot(luca, filippo, pch=19, col="blue", cex=2, xlab="rubbish", ylab="biomass")

# Installing packages
# CRAN
install.packages("terra")
library(terra)

install.packages("devtools")
library(devtools)

install_github("ducciorocchini/imageRy")
library(imageRy)

# Rtools 4.4
# update packages: avoid! just press enter


# installato imageRy studia bene le funzioni che voglio mettere nel progetto
# anche terra bisogna studiarselo bene perchè è enorme e Rocio ha detto che è super consigliato e utile
# riguarda anche vecchi appunti dell'esame di R e statistica che ci sono tutti i codici di ggplot2 per i grafici e anche altro
# salvare i pdf con dev.off() MEMO






