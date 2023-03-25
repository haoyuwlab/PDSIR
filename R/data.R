#' Weekly infection counts in Gueckedou
#'
#' Weekly infection counts of Ebola in the prefecture Gueckedou in Guinea during the 2013-2025 outbreak.
#' The full data for Guinea can be found in `ebola_guinea`.
#'
#' @format ## `ebola_gueckedou`
#' A list of three vectors:
#' \describe{
#'   \item{I_k}{number of weekly infections}
#'   \item{dates}{date of the week}
#'   \item{ts}{endpoints of observation intervals (units: days), used in the MCMC algorithm; the intervals correspond to the weeks in `dates`}
#' }
"ebola_gueckedou"

#' Weekly infection counts in each prefecture of Guinea
#'
#' These are the raw data for Guinea
#'
#' @format ## `ebola_guinea_raw`
#' A list of three vectors:
#' \describe{
#'   \item{I_k}{number of weekly infections}
#'   \item{dates}{date of the week}
#'   \item{ts}{endpoints of observation intervals (units: days), used in the MCMC algorithm; the intervals correspond to the weeks in `dates`}
#' }
"ebola_guinea_raw"
