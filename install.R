# Install the development version from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("wenbostar/PGA")
install.packages("tidyverse")
install.packages("optparse")
