#  Generic install script for proximity labeling pipeline
install.packages("reshape")
install.packages("optparse")
install.packages("stringr")
install.packages("ggplot2")
install.packages("extrafont")
install.packages("ggrepel")
install.packages("VennDiagram")
install.packages("data.table")
install.packages("readr")

library("extrafont")
loadfonts()

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("limma")