
<!-- README.md is generated from README.Rmd. Please edit that file -->

# retel-paper

Reproducible code for generating figures and outputs from simulations
and the application presented in the paper ‘Regularized Exponentially
Tilted Empirical Likelihood for Bayesian Inference,’ accessible on
[arXiv](https://arxiv.org/abs/2312.17015).

The following sections outline the directory structure and provide
instructions for replicating the results.

## Directory Structure

    .
    ├── application/        <- 📂 Application outputs
    ├── code/               <- 📂 All code for reproducibility
    │   ├── application/        <- 📂 Code for application in Section 6.1
    │   ├── figures/            <- 📂 Code for figures 
    │   ├── simulations/        <- 📂 Code for simulations
    ├── renv/               <- 📂 R environment and package dependency management
    ├── simulations/        <- 📂 Simulation outputs
    ├── .Rprofile           <- 🔧 Configuration file for R environment, including settings for renv
    ├── README.Rmd          <- 📝 README file written in R Markdown format
    ├── README.md           <- 📝 README file generated from README.Rmd
    ├── renv.lock           <- 📦 Lock file to set up the environment
    ├── retel-paper.Rproj   <- 🖥️ R Project file for RStudio 

## Instructions for Reproducibility

Please ensure Git, R (version 4.1.0 or higher), and preferably RStudio
(for interactive script execution on a personal computer) are installed.

1.  Clone this repository.

2.  Set the working directory to the `retel-paper/` folder or open the
    `retel-paper.Rproj` project file in RStudio.

3.  Restore the project library by executing the following command in R:

``` r
# install.packages("renv")
renv::restore()
```

This command will automatically install all the required packages,
including ‘retel.’ To maintain reproducibility and prevent potential
issues, it is strongly advised not to manually install the packages.

4.  Run the R scripts in the `code/` folder. Each script is designed to
    run independently, with script names corresponding to specific
    methods, configuration settings, or outputs from the paper. No edits
    or modifications are needed.

For example:

- `code/application/retel_f.R` reproduces the RETEL_f results for Table
  4.
- `code/figures/figure1.R` reproduces Figure 1.
- `code/simulations/cr/l0/n5s05.R` reproduces the results for Table 2
  when n = 5 and s = 0.5.

Note: The simulations and application code were run in a
high-performance computing (HPC) environment. The total run time for
sequentially executing all scripts in the `code/` folder can vary,
potentially extending to many hours, depending on the machine or
hardware specifications.
