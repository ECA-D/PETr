#' Vapour pressure
#' @description This function computes vapour pressure from temperature
#'
#' @param tm  = Temperature    (degrees Celcius)
#' @param vp  = Vapour pressure    (hPa)
#' @return A vector containing the daily vp
#'
vapour_pressure <- function(tm){
  vp = (17.38*tm)/(239.0+tm)
  vp = 6.107*exp(vp)
  return(vp)
}
