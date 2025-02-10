# Bioinformatics TCGA Breast Cancer Analysis

## Overview
This project focuses on analyzing breast cancer data from The Cancer Genome Atlas (TCGA) using R and Bioconductor. The analysis aims to identify key biomarkers for early detection and prognosis prediction.

## Features
- Gene expression analysis using Bioconductor packages (e.g., DESeq2, edgeR)
- Integration of TCGA data using `TCGAbiolinks`
- Principal Component Analysis (PCA) and visualization
- Machine learning models for cancer subtype classification

## Folder Structure
```
BIN_TCGA_BCA/
│-- TCGA_results/         # Processed results and outputs
│-- TCGA_data/            # Raw TCGA data (not included in the repo)
│-- megatron.R            # Main R script for analysis
│-- README.md             # Project documentation
│-- .gitignore            # Ignoring large files like TCGA_data/
```

## Installation & Setup
### 1. Clone the Repository
```bash
git clone https://github.com/tahagill/Bioinformatics_TCGA_Breast_Cancer_Analysis.git
cd Bioinformatics_TCGA_Breast_Cancer_Analysis
```

### 2. Install Required Packages
Open R and run the following:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "DESeq2", "edgeR", "ggplot2"))
```

### 3. Download TCGA Data
Since the TCGA dataset is large, it is not included in this repository. You can download it using the script in `megatron.R`. Open R and run:
```r
source("megatron.R")
```
This will download the necessary TCGA data and store it in `TCGA_data/`.

## Running the Analysis
Run the R script from the terminal or within RStudio:
```r
source("megatron.R")
```
The processed results will be saved in the `TCGA_results/` folder.

## Contributing
Feel free to fork the repository and submit pull requests!


