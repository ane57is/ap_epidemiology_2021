# Antonine Plague

This repository contains the code and documentation for a computational historiographical study of the Antonine Plague. The study employs Ordinary Differential Equation (ODE) models to explore various hypothesized transmission routes of disease spread and conducts sensitivity analyses to assess the impact of different parameters.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13857741.svg)](https://doi.org/10.5281/zenodo.13857741)


## Repository Contents

### R Scripts and Functions

- `Plague_model_functions.R`: This script sets up a system of ODEs for SIR and SEIR models of plague, smallpox, and measles, using the `deSolve` library in R.
- `PlottingTimeCourses.R`: This script plots results from the SIR and SEIR models, using `deSolve`, `tidyr`, and `ggplot2` libraries.

### Sensitivity Analysis

- `LHSnonuniform.R`: This script creates uniform and non-uniform LHS distributions for PRCC analysis of SIR and SEIR models.
- `GlobalSensitivityAnalysis.R`: This script conducts global sensitivity analysis for each model.

### R Markdown Files

- `AntoninePlagueModelingFigures.Rmd`: This R Markdown file produces time courses of ODE models defined in `Plague_model_functions.R` and generates Figure 1 as `Fig1_all_SEIR.tiff`, and Supplementary Figures S1 as `FigS1_plg_all_SEIR.tiff` and S2 as `FigS2_plg_all_SEIR.tiff`.
- `UniformLHSPRCC.Rmd`: Using uniform LHS distributions, this file performs sensitivity analysis and creates Supplementary Figure S3 as `FigS3.tiff`.
- `NonUniformLHSPRCC.Rmd`: Using non-uniform LHS distributions, this file performs sensitivity analysis with non-uniform parameter distributions and produces Supplementary Figures S4 to S8 as files `FigS4.tiff`, `FigS5.tiff`, `FigS6.tiff`, `FigS7.tiff`, and `FigS8.tiff`.
