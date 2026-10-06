# README file

The scripts in this repository can be used to reproduce the simulation results in *Spatio-temporal quantile regression forests for vegetation productivity in Northern Italy* by L. Merlo, L. Petrella, B. Foroni, V. L. Sciabolazza, and L. Salvati.

## Prerequisites
### Software requirements

-   [R](https://cran.r-project.org/) version 4.5.3 or higher
-   [RStudio](https://rstudio.com/) version 2026.07.1+147 or higher

The code has been parallelized using the following packages:
-   `foreach`
-   `parallel`
-   `doParallel`

## Script description
-    `mixed_model_function_time_space.R` contains the main functions to implement the proposed method.
-    `PIRLS_time_space.R` contains the iterative reweighted least squares algorithm.


## Simulation study
Run the `time_space_server.R` scripts to reproduce the results in Section S1 of the Supplementary Materials (Tables S1 and S2).
The scenario can be set by changing the object `scenario` to "linear" or "non_linear".
