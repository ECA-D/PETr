#' Potential Evapotranspiration
#' @description This function takes a dataframe object as input and computes PET.
#'
#' @param indat = dataframe containing the input variables
#'
#' @return A vector containing the daily Potential EvapoTranspiration (mm/day)
#' @examples
#' pet_penmon(indat)
#'
#' Where the dataframe "indat" contains the following variables:
#'
#' tmin  = Temperature min    (degrees Celsius)
#' tmax  = Temperature max    (degrees Celsius)
#' only one of the following variables (vp OR rh OR dp) are to be passed
#' vp  = Vapour pressure    (hPa)
#' rh =  Relative humidity  (%)
#' dp =  Dewpoint           (degrees Celsius)
#' ws  = Wind speed at 10m  (m/s)
#' either ss or rs is passed, but not both
#' ss  = Sunshine duration  (0-24 hours)
#' cl  = Cloudiness         (0-1)
#' cl is optional for use with ss, if unavailble cl is estimated
#' rs  = radiation          (kwh/m^2)
#' Also rquired is the latitude and longitude
#' lat = Latitude  (degrees)
#' lons = Longitude (degrees)
#' Elevation (needed) is retrieved from the PETr internal file "data/elev_dat.rda" which contains data from:
#' http://www.ecad.eu/download/ensembles/data/Grid_0.1deg_reg_ensemble/elev_ens_0.1deg_reg_v17.0e.nc
#' Missing values should be converted to NA
#'
#' @export

pet_penmon <- function(indat) {
  data("elev_dat", package="PETr")
#browser()
  lat <- indat$lat
  lon <- indat$lons
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

  jd <- as.POSIXlt(indat$tmax.dates)$yday + 1
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
  }

  # check for pressence of ss and cl
  if (is.null(indat$ss)) {
    ss_pres <- FALSE
  } else {
    ss_pres <- TRUE
  }

  if (is.null(indat$cl)) {
    cld_pres <- FALSE
  } else {
    cld_pres <- TRUE
  }

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
  psy <- 0.067
  # estimate 2m ws
  wnf <- log((2.0 - 0.08) / 0.015)
  wnf <- wnf / log((10.0 - 0.08) / 0.015)
  wn2m <- wnf * indat$ws
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

  cl <- array(dim=length(indat$tmax))
  if(ss_pres) {

    if (!cld_pres) {
      # estimate cloudiness (cl)
      dyln <- 7.64 * snst
      ss_gt <- which(indat$ss > dyln)
      indat$ss[ss_gt] <- dyln[ss_gt]
      dyln_zero <- which(dyln < 0.000001)
      dyln_nonzero <- which(dyln >= 0.000001)
      cl[dyln_zero] <- 1.0
      cl[dyln_nonzero] <- 1.0 - indat$ss[dyln_nonzero] / dyln[dyln_nonzero]
    }
    ssf <- 1.0 - cl

    rns <- 0.77 * (0.25 + (0.50 * (ssf))) * dter
    rnl <- sigma * (0.9 * (ssf) + 0.1) * (0.34-0.14*sqrt(avp)) * (ktmin^4.0 + ktmax^4.0)
  } else {

    rso <- (0.75 + (2 * 10^-5) * elev) * dter
    rns <- (1 - a) * indat$rs
    rs_div_rso <- indat$rs / rso
    rs_div_rso[rs_div_rso > 1.0] <- 1.0
    rnl <- sigma * (0.34 - 0.14 * sqrt(avp)) * (ktmax^4 + ktmin^4)/2 * (1.35 * rs_div_rso - 0.35)
    rnl[rnl < 0] <- 0
  }
  # Net radiation
  rnet <- rns - rnl
  vptmin <- 0.1 * vapour_pressure(indat$tmin)
  vptmax <- 0.1 * vapour_pressure(indat$tmax)
  vpd <- ((vptmax + vptmin) / 2.0) - avp
#browser()
  pet <- 0.408 * dvp * (rnet - tflux)
  pet <- pet + (900.0 / (tm + 273.16)) * (psy * wn2m * vpd)
  pet <- pet / (dvp + (psy * (1.0 + (wn2m * 0.34))))

  pet[pet < 0] <- 0
  return(pet)
}
#' Vapour pressure
#' @description This function computes vapour pressure from temperature
#'
#' @param tm  = Temperature    (degrees Celcious)
#' @param vp  = Vapour pressure    (hPa)
#' @return A vector containing the daily vp
#'
vapour_pressure <- function(tm, vp){
  vp = (17.38*tm)/(239.0+tm)
  vp = 6.107*exp(vp)
  vapour_pressure = vp
}
#' Relative Humidity
#' @description This function computes RH from dewpoint
#'
#' @param t  = Temperature mean    (degrees Celcious)
#' @param dp  = Dewpoint    (degrees celsius)
#' @return A vector containing the daily vp
#'

dewpoint_to_humidity <- function (dp = NA, t = NA)
{
  if (length(dp) != length(t)) {
    stop("The vectors for temperature('t') and dewpoint temperature ('dp') must have the same length.")
  }
  if (length(dp[dp > t & !is.na(dp) & !is.na(t)]) > 0) {
 #   browser()
    dp[dp > t] <- t[dp > t]
  }
  beta <- (112 - (0.1 * t) + dp) / (112 + (0.9 * t))
  relative.humidity <- beta^8
  return(relative.humidity)
}
