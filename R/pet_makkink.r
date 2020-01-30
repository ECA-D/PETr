#' Potential Evapotranspiration
#' @description This function takes a dataframe object as input and computes PET using Makkink's method.
#'
#' @param indat = dataframe containing the input variables
#'
#' @return A vector containing the daily Potential EvapoTranspiration (mm/day)
#' @examples
#' pet_makkink(indat)
#'
#' Where the dataframe "indat" contains the following variables:
#'
#' tmin  = Temperature min    (degrees Celsius)
#' tmax  = Temperature max    (degrees Celsius)
#' rs  = radiation          (MJ/m^2)
#'
#' Missing values should be converted to NA
#'
#' @export

pet_makkink <- function(indat) {
  data("elev_dat", package="PETr")
  tm <- rowMeans(cbind(indat$tmin, indat$tmax))
  ndat <- length(tm)

  # line 26-37 is just to get land only
  lat <- indat$lat
  lon <- indat$lon
  # lat index
  lat_index <- which.min(abs(lat - elev_lat))
  # lon index
  lon_index <- which.min(abs(lon - elev_lon))
  elev <- elev_dat[lon_index, lat_index]
  if(is.na(elev))
  {
    pet <- array(dim=ndat, NA)
    return(pet)
  }
  Q <- indat$rs
  vptm <- 0.1 * vapour_pressure(tm)
  dvp <- (4099.0 * vptm) / (tm + 237.3)^2.0
  # psy calc
  lambda <- 2.45
  p <- 101.3*((293-(0.0065*elev))/293)^5.26
  psy <- 0.00163*(p/lambda)
  dpsy <- dvp / (dvp+psy)
  pet_mk <- 0.61 * (dpsy * Q/2.45) - 0.12

  return(pet_mk)
}
