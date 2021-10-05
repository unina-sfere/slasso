## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib slasso, .registration = TRUE
## usethis namespace: end
NULL



#' @title slasso
#' @details
#'
#'\tabular{ll}{
#'Package: \tab rofanova\cr
#'Type: \tab Package\cr
#'Version: \tab `r packageVersion("rofanova")` \cr
#'Date: \tab  `r Sys.Date()` \cr
#'License: \tab `r packageDescription("rofanova", fields="License")`\cr
#'}
#'
#'
#'
#'
#' @aliases {rofanova}-package
#' @author Fabio Centofanti, Antonio Lepore, Biagio Palumbo
#' @references
#' Centofanti, F., Lepore, A., & Palumbo, B. (2021).
#'
#' @seealso \code{\link{fusem}} \code{\link{funmad}}
#' @examples
##' library(rofanova)
##' library(rofanova)
#' data_out<-simulate_data(scenario="one-way")
#' label_1=data_out$label_1
#' X_fdata<-data_out$X_fdata
#' B=10
#' cores=1
#' per_list_median<-rofanova(X_fdata,label_1,B = B,family="median",cores=cores)
#' pvalue_median_vec<-per_list_median$pval_vec
#' per_list_huber<-rofanova(X_fdata,label_1,B = B,family="huber",cores=cores)
#' pvalue_huber_vec<-per_list_huber$pval_vec
#' per_list_bisquare<-rofanova(X_fdata,label_1,B = B,family="bisquare",cores=cores)
#' pvalue_bisquare_vec<-per_list_bisquare$pval_vec
#' per_list_hampel<-rofanova(X_fdata,label_1,B = B,family="hampel",cores=cores)
#' pvalue_hampel_vec<-per_list_hampel$pval_vec
#' per_list_optimal<-rofanova(X_fdata,label_1,B = B,family="optimal",cores=cores)
#' pvalue_optimal<-per_list_optimal$pval
#'@import fda.usc robustbase
"_PACKAGE"



