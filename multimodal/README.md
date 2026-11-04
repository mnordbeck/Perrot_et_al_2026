# Multiblock analysis of clinical and biological measurements in a mouse model of multiple sclerosis

This repository contains the R scripts used to preprocess the data,
perform multiple co-inertia analysis (MCOA), and investigate associations
between MCOA scores and donor characteristics.

## Data structure

Each row represents one mouse. Multiple mice may receive cells from the
same human donor, identified by `Subject_ID`.

The analyses combine biological measurements and clinical measurements
collected at the 3-month time point. Donor metadata include donor group
(HD: healthy donor; MS: multiple sclerosis donor), age, sex, and treatment.

## Files

- `main.R`: Runs the analysis pipeline.
- `generate_data.R`: Reads and preprocesses the source data, retaining
  the biological variables and clinical measurements at 3 months.
  This script is called by `rgcca_analysis.R`.
- `rgcca_analysis.R`: Performs MCOA using the RGCCA package and displays:
  - bootstrap results for the first component;
  - a biplot;
  - donor barycenters, calculated by averaging mouse scores within donors;
  - donor barycenter plots colored by age, sex, and treatment.
  It also creates the `dat` object used by `effect_analysis.R`.
- `effect_analysis.R`: Fits linear mixed-effects models relating the
  first MCOA component score to donor treatment, sex, or age.

## Input data

- Source dataset: `Subject_Tables.xlsx`



## Software requirements

The analysis requires R and the following packages (tested with R 4.5.0):

- `readxl`: importing Excel data. (tested with version 1.5.0.1)
- `RGCCA`: multiple co-inertia analysis and bootstrap analysis. (tested with version 3.0.3)
- `ggplot2`: visualization. (tested with version 4.0.3)
- `lmerTest` : linear mixed-effects models and statistical tests. (tested with version 3.2.1)
- `dplyr`: data manipulation. (tested with version 1.2.1)
- `grid`: graphical elements, including arrows and units. (included in base R)

Install the required packages once:

    install.packages(c("readxl", "dplyr", "RGCCA", "ggplot2", "lmerTest"), dependencies = TRUE)

## Running the analysis

Set the working directory to the directory containing `main.R`
and run:

    source("main.R")

Alternatively, from a terminal in that directory:

    Rscript main.R

The execution order is:

1. `rgcca_analysis.R`, which calls `generate_data.R`.
2. `effect_analysis.R`, which uses the objects created in the first step.

## Statistical analyses

MCOA is performed using `RGCCA::rgcca()` with `method = "mcoa"`
and two components.

The mixed-effects analyses are restricted to MS donors with recorded
treatment information. Three separate models are fitted:

    comp1 ~ Treatment + (1 | Subject_ID)
    comp1 ~ Sex + (1 | Subject_ID)
    comp1 ~ Age_c + (1 | Subject_ID)

`Age_c` is donor age centered at 40 years. The treatment reference
category is `"None"` (no treatment), and the sex reference category
is `"F"` (female).

Each model includes a random intercept for donor to account for
dependence among mice sharing the same donor. The three models are
not mutually adjusted: each includes only one donor characteristic
as a fixed effect.

Fixed-effect tests use the Satterthwaite approximation implemented
in `lmerTest`.

## Outputs

Processed data and analysis results remain available as R objects
in the current session after running `source("main.R")`.

Model summaries and bootstrap results are printed to the console.
Figures are displayed in the active R graphics device.

No explicit export of tables or figures is performed by the scripts.



## References

- Girka F, Camenen E, Peltier C, Gloaguen A, Guillemot V,
  Le Brusquet L, Tenenhaus A (2023).
  *RGCCA: Regularized and Sparse Generalized Canonical Correlation
  Analysis for Multiblock Data*. R package version 3.0.3.
  https://CRAN.R-project.org/package=RGCCA


- Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017).
  *lmerTest Package: Tests in Linear Mixed Effects Models*.
  Journal of Statistical Software, 82(13), 1–26.
  https://doi.org/10.18637/jss.v082.i13



