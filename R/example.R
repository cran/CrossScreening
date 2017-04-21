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

#' Lead in children
#'
#' Morton et al. (1982) compared the levelsof lead in the blo od of 33 children wh ose fathers worked in afactory w here lead was used in the manufacture of batteriesto 33 lead levels in ma tched control ch ild ren of the sam e agefrom the same neighborhood. The variables are as follows:
#'
#' \describe{
#' \item{control}{lead levels (ug/dl)}
#' \item{level}{father's potential exposure}
#' \item{hyg}{hygine of father employed in the lead industry}
#' }
#'
#' @references Morton, D. E., Saah, A. J., Silberg, S. L., Owens, W. L., ROBERTS, M. A., & Saah, M. D. (1982). Lead absorption in children of employees in a lead-related industry. American Journal of Epidemiology, 115(4), 549-555.
#'
#' @docType data
#'
#' @usage data(lead)
#'
#' @format A data.frame.
#'
#' @keywords datasets
#'
"lead"

#' Methotrexate workers
#'
#' Methotrexate is a drug used totreat cancer, but there is concern that it may itself be carcinogenic in healthy individuals who are exposed while either manufacturing or administering the drug.
#' Deng et al. (2005) compared 21 workers engaged in the production of methotrexate to 21 unexp osed controls matched for age, gender, and smoking. The variables are (prefix "w" means exposed and "c" means control)
#' \describe{
#' \item{Mftcr}{mutant frequency of TCR gene}
#' \item{Mfhrpt}{mutant frequency of hprt gene}
#' \item{mtl}{mean tail length}
#' \item{mtm}{mean tail moment}
#' \item{id}{identifier}
#' \item{sex}{sex}
#' \item{age}{age}
#' \item{smoke}{smoking}
#' \item{years}{exposure years}
#' \item{protection}{protection measures, G for gloves, M for mask, N for none}
#' \item{mask}{if the protection includes mask}
#' }
#'
#' @docType data
#'
#' @usage data(methotrexate)
#'
#' @format A data.frame.
#'
#' @keywords datasets
#'
#' @references Deng, H., Zhang, M., He, J., Wu, W., Jin, L., Zheng, W., ... & Wang, B. (2005). Investigating genetic damage in workers occupationally exposed to methotrexate using three genetic end-points. Mutagenesis, 20(5), 351-357.
#'
"methotrexate"

#' Obtains treatment-minus-control differences in the \code{nhanes.fish} dataset
#'
#' @export
#'
#' @return a 234 * 46 matrix of the log2 differences
#'
nhanes.log2diff <- function() {
    data("nhanes.fish", package = "CrossScreening", envir=environment())
    data("nhanes.fish.match", package = "CrossScreening", envir=environment())
    data <- get("nhanes.fish", envir=environment())
    match <- get("nhanes.fish.match", envir=environment())
    outcomes <- grep("^o\\.", names(data))
    log2diff <- function(y1, y2) {
        if (min(c(y1, y2)) == 0) {
            y1 <- y1 + 1
            y2 <- y2 + 1
        }
        log2(y1) - log2(y2)
    }
    d <- data.frame(sapply(outcomes, function(j) log2diff(data[match$treated, j], data[match$control, j])))
    colnames(d) <- names(data)[outcomes]
    d
}
