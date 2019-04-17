##Post-pipeline RNA-seq

###Install/update R and RStudio:

- Update R to 3.5.3: https://cran.r-project.org/ (use link at top for *precompiled binary*. **Skip if you already have R version 3.5.something**
- Install RStudio (if you don't have it already installed): https://www.rstudio.com/products/rstudio/download/. **NOTE: install the FREE version.**

###R packages to install:

Open RStudio and install the packages we'll use by running these commands below:

- magrittr: `install.packages("magrittr", dependencies = TRUE)`
- statmod: `install.packages("statmod", dependencies = TRUE)`
- BiocManager: `install.packages("BiocManager")`
- edgeR: `BiocManager::install("edgeR", version = "3.8")`
- org.Mm.eg.db: `BiocManager::install("org.Mm.eg.db", version = "3.8")`
