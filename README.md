# Feature Selection and FDR Control in High-Dimensional Data

This repository contains a collection of R scripts designed for feature selection and False Discovery Rate (FDR) control in high-dimensional data. The scripts implement various methods such as PMDC-DS, PMDC-MDS, REDS, and Knockoff-PMDC, and are used for simulations and analyses in drug resistance studies.

## Table of Contents
1. [Overview](#overview)
2. [Script Descriptions](#script-descriptions)
3. [Usage](#usage)
4. [Dependencies](#dependencies)
5. [Data Files](#data-files)
6. [License](#license)

## Overview
The scripts in this repository are designed to perform feature selection and FDR control in high-dimensional datasets. They are particularly useful in scenarios where the number of features (e.g., genetic markers) is much larger than the number of samples. The methods implemented include:

- **PMDC-DS**: Projection-based Martingale Difference with Data Splitting.
- **PMDC-MDS**: Projection-based Martingale Difference with Multiple Data Splitting.
- **REDS**: Robust Estimation with Data Splitting.
- **Knockoff-PMDC**: Knockoff framework combined with PMDC for FDR control.
- **PMDD-BHq**: Application of PMDC statistics to the Benjamini-Hochberg procedure.

## Script Descriptions
The scripts are categorized based on their functionality:

### 1. **Main Scripts**
- **Case 1**
  - **`HIV_drug.R`**: Processes drug resistance data, performs feature selection using the Knockoff method, and calculates FDR.
- **Case 2**
  - **`GSE5680_dsscreen.R`**: Feature screening using PMDC-DS and REDS.
  - **`GSE5680_mdsscreen.R`**: Feature screening using PMDC-MDS.
  - **`GSE5680_BHq.R`**: Feature screening using PMDD-BHq.
- **`main.R`** (RUN ME): Main simulation script that sets up parameters and runs simulations for different models and methods via data splitting.
- **`MDC_HT.R`**: Runs simulations for the MDC method with different hard thresholding techniques.

### 2. **Component Files of Feature Screening Methods**
- **`DSscreen.R`**: Implements the PMDC-DS (Data Splitting) screening procedure.
- **`MDSscreen.R`**: Implements the PMDC-MDS (Multiple Data Splitting) screening procedure.
- **`DSBHq.R`**: Implements the PMDD-BH (Benjamini-Hochberg) screening procedure.
- **`REDScreen.R`**: Implements the REDS (REflection via Data Splitting) screening procedure.
- **`mirrorstat_pmdc.R`**: Computes mirror statistics for PMDC-DS methods.
- **`knockoff_pmdc.R`**: Computes mirror statistics for the Knockoff-PMDC method.
- **`REDS_pmdc.R`**: Computes mirror statistics for the REDS-PMDC method.
- **`pmdd_H.R`**: Computes the mean and variance of PMDD.

### 3. **Supporting Files**
- **Main Support Files**
  - **`generate_data.R`**: Generates synthetic data for simulations based on different models (e.g., linear, exp-additive).
  - **`PMDC-DS.R`**: Runs the simulation for the PMDC-DS method.
  - **`PMDC-MDS.R`**: Runs the simulation for the PMDC-MDS method.
  - **`PMDD-BHq.R`**: Runs the simulation for the PMDD-BH method.
  - **`REDS.R`**: Runs the simulation for the REDS method.
  - **`Knockoff-PMDC.R`**: Runs the simulation for the Knockoff-PMDC method.
  - **`Screen_FDR.R`**: Aggregates screening functions.
  - **`result_print.R`**: Prints the results of simulations.
- **Case 1**
  - **`drug_result.R`**: Prints the results of the HIV drug resistance case.
- **Case 2**
  - **`GSE5680_screen.R`**: Selects the top-ranked genes as the ground truth.

## Usage
To run the simulations or analyses, follow these steps:

1. **Set up the environment**: Ensure you have R installed along with the necessary libraries (`mvtnorm`, `MASS`, `knockoff`, `Matrix`).
2. **Run the main script**: Execute `main.R` to start the simulation. Modify the parameters in `main.R` to customize the simulation settings.
   
   **Example command:**
   ```r
   source("main.R")
   ```
   
3. **Analyze the results**: The results are saved as `.RData` files, which can be loaded and analyzed in R.

## Dependencies
- R (>= 3.6.0)
- Libraries: `mvtnorm`, `MASS`, `knockoff`, `Matrix`

## Data Files
- **Case 1**
  - **`TSM.csv`**: Ground truth for pinpointing factors related to drug resistance.
  - **`PI_DATA.csv`**: Covariates corresponding to PI treatment.
  - **`NRTI_DATA.csv`**: Covariates corresponding to NRTI treatment.
  - **`phs000276_screen`**: Categorized genotype data according to key pinpoints.
- **Case 2**
  - **`GSE5680.Rdata`**: Covariates for Case 2 (X).
  - **`18976_x.Rdata`**: Response variable for Case 2 (Y).
  - **`GSE5680_05.Rdata`**: Top-ranked active features for τ=0.5 (i.e., stored results from `GSE5680_screen.R`).

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

For any questions or issues, please open an issue in the repository or contact the maintainers.

