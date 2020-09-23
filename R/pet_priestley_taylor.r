#' Potential Evapotranspiration
#' @description This function takes a dataframe object as input and computes PET using Makkink's method.
#'
#' @param indat = dataframe containing the input variables
#'
#' @return A vector containing the daily Potential EvapoTranspiration (mm/day)
#' @examples
#' pet_priestley_taylor(indat)
#'
#' Where the dataframe "indat" contains the following variables:
#'
#' tmin  = Temperature min    (degrees Celsius)
#' tmax  = Temperature max    (degrees Celsius)
#' only one of the following variables (vp OR rh OR dp) are to be passed
#' vp  = Vapour pressure    (hPa)
#' rh =  Relative humidity  (%)
#' dp =  Dewpoint           (degrees Celsius)
#' rs  = radiation          (MJ/m^2)
#' Also rquired is the latitude and longitude (for land sea mask)
#' lat = Latitude  (degrees)
#' lons = Longitude (degrees)
#'
#' Missing values should be converted to NA
#'
#' @export

pet_priestley_taylor <- function(indat) {
  data("elev_dat", package="PETr")
#  browser()
  lat <- indat$lat
  lon <- indat$lon
  # lat index
  lat_index <- which.min(abs(lat - elev_lat))
  # lon index
  lon_index <- which.min(abs(lon - elev_lon))
  tm <- rowMeans(cbind(indat$tmin, indat$tmax))
  ndat <- length(tm)

  elev <- elev_dat[lon_index, lat_index]
  if(is.na(elev))
  {
    pet <- array(dim=ndat, NA)
    return(pet)
  }
  # check whether vp, rh or dp are present
  if (is.null(indat$vp) && !is.null(indat$rh)) {
    # estimate vp from rh
    vps <- vapour_pressure(tm)
    vp <- (indat$rh / 100) * vps

  } else if (is.null(indat$vp) &&  !is.null(indat$dp)) {
    # estimate vp from dp
    vps <- vapour_pressure(tm)
    rh <- dewpoint_to_humidity(indat$dp, tm)
    vp <- (rh / 100) * vps

  }  else if (!is.null(vp)) {
    vp <- indat$vp
  }else {
    stop("Terminating - missing imput variables")
  }
  # calc day count
  if(is.null(indat$jd))
  {
    jd <- as.POSIXlt(indat$tmax.dates)$yday + 1
  } else {
    jd <- indat$jd
  }
  # real ssf      # Sunshine factor
  # real slrd     # Solar declination
  # real snst     # Sunset hour angle
  # real cos_snst # Cosinus of (Sunset hour angle)
  # real dr       # Relative distance Earth to Sun
  # real dter     # Daily total extraterrestrial radiation
  # real dvp      # Slope of the vapour pressure
  # real tflux    # Daily temperature fluctuations
  # real rnet     # Net radiation
  # real dyln     # Day length
  # real psy      # Psychrometric constant set to the same as in SWURVE 0.067 kPaC-1

  ptm <- array(dim=ndat)
  ptm[1] <- tm[1]
  ptm[2:ndat] <- tm[1:(ndat-1)]
  a <- 0.23 # albedo
  sigma <- 2.45e-09
  pi <- 3.1415
  avp <- 0.1 * vp
  vptm <- 0.1 * vapour_pressure(tm)
  dvp <- (4099.0 * vptm) / (tm + 237.3)^2.0
  tflux <- 0.38 * (tm - ptm)
  # psy calc
  lambda <- 2.45
  p <- 101.3*((293-(0.0065*elev))/293)^5.26
  psy <- 0.00163*(p/lambda)
  # convert lat to radians
  latr <- lat * (pi / 180.0)
  # Solar declination
  slrd <- 0.409 * sin((0.0172 * jd) - 1.39)
  # Relative distance Earth to Sun
  dr <- 1.0 + 0.033 * cos(0.0172 * jd)
  # solar constant
  Gsc <- 0.0820
  ktmin <- indat$tmin + 273.16
  ktmax <- indat$tmax + 273.16
  # deals with S. Hemisphere...
  cos_snst <- (-1) * tan(latr) * tan(slrd)
  cos_snst[cos_snst > 1] <- 1
  cos_snst[cos_snst < -1] <- -1
  # Sunset hour angle
  snst <- acos(cos_snst)
  # Daily total extraterrestrial radiation
  dter <- 37.6 * dr * (snst * sin(latr) * sin(slrd) + cos(latr) * cos(slrd) * sin(snst))

  Q <- indat$rs
  rso <- (0.75 + (2 * 10^-5) * elev) * dter
  rns <- (1 - a) * Q
  rs_div_rso <- Q / rso
  rs_div_rso[rs_div_rso > 1.0] <- 1.0
  rnl <- sigma * (0.34 - 0.14 * sqrt(avp)) * (ktmax^4 + ktmin^4) * (1.35 * rs_div_rso - 0.35)
  rnl[rnl < 0] <- 0
  dpsy <- dvp / (dvp+psy)
  # net incoming shortwave radiation
  rns <- (1 - a) * Q
  # Net radiation
  rnet <- rns - rnl

#  C <- 1.26  # original constant
  C <- 1.115
  pet_pt = C*dpsy*(rnet-tflux)/2.45
  pet_pt[pet_pt < 0] <- 0
  return(pet_pt)
}
