
<!-- README.md is generated from README.Rmd. Please edit that file -->

# retel-paper

Reproducible code for generating figures and outputs from simulations
and the application presented in the paper ‘Regularized Exponentially
Tilted Empirical Likelihood for Bayesian Inference,’ accessible on
[arXiv](https://arxiv.org/abs/2312.17015). The following sections
outline the directory structure and provide detailed instructions for
replicating the results.

## Directory Structure

    .
    ├── application/        <- 📂 Application outputs
    ├── code/               <- 📂 All code for reproducibility
    │   ├── application/        <- 📂 Code for application in Section 5
    │   ├── figures/            <- 📂 Code for figures 
    │   │   ├── supplement/         <- 📂 Supplement figures  
    │   ├── simulations/        <- 📂 Code for simulations
    │   │   ├── cr/                 <- 📂 Simulations in Section 4.1
    │   │   │   ├── l0/                 <- 📂 Table 2 when l = 0
    │   │   │   ├── l2/                 <- 📂 Table 3 when l = 2
    │   │   ├── kl/                 <- 📂 Simulations in Section 4.2
    │   │   │   ├── el/                 <- 📂 EL implementation
    │   │   │   ├── etel/               <- 📂 ETEL implementation
    │   │   │   ├── retel_f/            <- 📂 RETEL_f implementation
    │   │   │   ├── retel_r/            <- 📂 RETEL_r implementation
    │   │   ├── mb/                 <- 📂 Simulations in Section 4.1
    │   │   │   ├── tau1/               <- 📂 Table 1 when tau_n = 1
    │   │   │   ├── taulogn/            <- 📂 Table 1 when tau_n = log(n)
    ├── figures/            <- 📂 Figure files (PDFs and PNGs)
    │   ├── supplement/         <- 📂 Supplement figures 
    ├── renv/               <- 📂 R environment and package dependency management
    ├── simulations/        <- 📂 Simulation outputs
    │   ├── cr/                 <- 📂 Outputs in Section 4.1
    │   │   ├── l0/                 <- 📂 Table 2 when l = 0
    │   │   ├── l2/                 <- 📂 Table 3 when l = 2
    │   ├── kl/                 <- 📂 Outputs in Section 4.2
    │   │   ├── el/                 <- 📂 EL outputs
    │   │   ├── etel/               <- 📂 ETEL outputs
    │   │   ├── retel_f/            <- 📂 RETEL_f outputs
    │   │   ├── retel_r/            <- 📂 RETEL_r outputs
    │   ├── mb/                 <- 📂 Outputs in Section 4.1
    │   │   ├── tau1/               <- 📂 Table 1 when tau_n = 1
    │   │   ├── taulogn/            <- 📂 Table 1 when tau_n = log(n)
    ├── .Rprofile           <- 🔧 Configuration file for R environment, including settings for renv
    ├── README.Rmd          <- 📝 README file written in R Markdown format
    ├── README.md           <- 📝 README file generated from README.Rmd
    ├── renv.lock           <- 📦 Lock file to set up the environment
    ├── retel-paper.Rproj   <- 🖥️ R Project file for RStudio 

## Instructions for Reproducibility
