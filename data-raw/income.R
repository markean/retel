## Code to produce the `income` data
library(dplyr)

if (file.exists("data-raw/income.csv")) {
  raw <- read.csv("data-raw/income.csv")
} else {
  raw <-
    read.csv("https://github.com/markean/retel/raw/main/data-raw/income.csv")
}

# Create package data
income <- raw |>
  # Create a new variable: adjusted median income
  mutate(ami = (pci_1989 / pci_1979) * mi_1979) |>
  # Standardize variables to have mean = 0 and standard deviation = 1
  mutate_at(c("mi_1989"), ~ (scale(., center = TRUE, scale = TRUE) |>
    # Needed to maintain the variable name
    as.vector())) |>
  mutate_at(c("mi_1979"), ~ (scale(., center = TRUE, scale = TRUE) |>
    # Needed to maintain the variable name
    as.vector())) |>
  mutate_at(c("ami"), ~ (scale(., center = TRUE, scale = TRUE) |>
    # Needed to maintain the variable name
    as.vector()))

usethis::use_data(income, overwrite = TRUE)
