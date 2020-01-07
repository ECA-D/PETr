#' Relative Humidity
#' @description This function computes RH from dewpoint
#'
#' @param t  = Temperature mean    (degrees celcius)
#' @param dp  = Dewpoint    (degrees celsius)
#' @return A vector containing the daily vp
#'

dewpoint_to_humidity <- function (dp = NA, t = NA)
{
  #browser()
  if (length(dp) != length(t)) {
    stop("The vectors for temperature('t') and dewpoint temperature ('dp') must have the same length.")
  }
  if (length(dp[dp > t & !is.na(dp) & !is.na(t)]) > 0) {
    dp[dp > t] <- t[dp > t]
  }
  # relative.humidity <- 100*((112 - (0.1 * t) + dp) / (112 + (0.9 * t)))^8
  b <- 17.625
  c <- 243.04
  rh <- 100*exp((c*b*(dp-t))/((c+t)*(c+t*dp)))
  rh[rh < 0] <- 0
  rh[rh > 100] <- 100
  return(rh)
}
