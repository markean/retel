#' Median Income for 4-Person Families in the USA
#'
#' A dataset of median income for 4-person families by state.
#'
#' @usage data("income")
#' @format A data frame with 51 rows and 6 columns:
#' \describe{
#'   \item{state}{States, including the District of Columbia.}
#'   \item{mi_1979}{Estimated median income for 4-person families in 1979 (standardized).}
#'   \item{mi_1989}{Estimated median income for 4-person families in 1989 (standardized).}
#'   \item{pci_1979}{Per capita income in 1979.}
#'   \item{pci_1989}{Per capita income in 1989.}
#'   \item{ami}{Census median income in 1979, adjusted for per capita income (standardized).}
#' }
#' @source \href{https://www.census.gov/data/tables/time-series/demo/income-poverty/4-person.html}{https://www.census.gov/data/tables/time-series/demo/income-poverty/4-person.html}
#' @examples
#' data("income")
#' income
"income"
