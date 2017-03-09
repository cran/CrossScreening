#' Health effects of fish
#'
#' Data from NHANES (2013-2014) with 1107 observations and 87 variables. Variables whose name start with "o." are lab measurements (such as blood mercury) that can be used as outcomes. The demographics and background variables include gender, age, income, indicator for missing income, race, education, indicator for smoked ever, number of cigararttes smoked in the last month. Fish intakes in the last month (in servings) are summed up in the "fish" variable, which is used to create the binary indicator "fish.level".
#'
#' @docType data
#'
#' @usage data(nhanes.fish)
#'
#' @format A data.frame.
#'
#' @keywords datasets
#'
#'
"nhanes.fish"

#' Pair matching result
#'
#' Each row is a matched pair, the first/second entry is the id of low/high fish intake in the \code{nhanes.fish} data frame.
#'
#' @docType data
#'
#' @usage data(nhanes.fish.match)
#'
#' @format A data.frame.
#'
#' @keywords datasets
#'
#'
"nhanes.fish.match"
