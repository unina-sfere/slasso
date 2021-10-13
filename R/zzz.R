.onLoad<-function(libname, pkgname){
  
if (!identical(objective@system, R.version$system)) 
  {
  # Rcpp::loadModule("RcppArmadillo", TRUE)
  # loadNamespace("Rcpp")
  # aa<-inline::getPlugin("Rcpp")
  # print(aa)
  # aa<-inline::getPlugin("RcppArmadillo")
  # print(aa)
  # print(.libPaths())
  # print(pkgname)
  # print(libname)
  # library.dynam("RcppArmadillo", pkgname, "C:/Program Files/R/R-4.0.5/library")
  # requireNamespace("RcppArmadillo", quietly = TRUE)
 # install.packages("RcppArmadillo",repos = "http://cran.us.r-project.org")
  objective <- cxxfunctionplus2(methods::signature(), body=objective.body,
                                includes=objective.include, plugin="RcppArmadillo",save.dso=TRUE)
  utils::assignInNamespace("objective", objective, ns = "slasso")

  gradient <- cxxfunctionplus2(methods::signature(), body=gradient.body,
                               includes=gradient.include, plugin="RcppArmadillo",save.dso=TRUE)
  utils::assignInNamespace("gradient", gradient, ns = "slasso")
  # useDynLib(gradient)
  # useDynLib(objective)
  # 
  }
}

# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Smooth LASSO Estimator for the Function-on-Function Linear Regression Model")
# }
# inline::inlineCxxPlugin
