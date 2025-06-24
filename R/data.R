#' Example test data for dose response curve
#'
#' A small example dataset containing dose-response measurements for a test drug.
#'
#' @format A data frame with 9 rows and 5 variables:
#' \describe{
#'   \item{drug}{Character, drug name (all "testdrug")}
#'   \item{nm}{Character, dose concentration in nanomolar (e.g., "1000", "500", etc.) and MUST also include "blank" as a value}
#'   \item{rep1}{Numeric, replicate 1 response value}
#'   \item{rep2}{Numeric, replicate 2 response value}
#'   \item{rep3}{Numeric, replicate 3 response value}
#' }
#' @source Internal example data for package demonstration
"testdata"
