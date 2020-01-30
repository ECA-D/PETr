library(gridclimind)
library(ncdf4)

in.dir <- "input/"
author.data <- list(ERA ="5")
out.dir <- "PET_output/"

############ PET Penmom-Moteith ############
### using sunshine hours (ss) derived from cloud cover ###
pet.file <- "PET_penman_ss_testfile.nc"
out.file <- sprintf("%s%s", out.dir, pet.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"tn_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"tx_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"dp_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"ss_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"ws_0.25deg_regular_1979-2018_sub.nc"))

create.pet.penman.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE)

### using radiation (rs) ###
pet.file <- "PET_penman_rs_testfile.nc"
out.file <- sprintf("%s%s", out.dir, pet.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"tn_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"tx_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"dp_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"rs_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"ws_0.25deg_regular_1979-2018_sub.nc"))

create.pet.penman.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE)

############ PET Makkink (MK) ############
pet.file <- "PET_makkink_testfile.nc"
out.file <- sprintf("%s%s", out.dir, pet.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"tn_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"tx_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"rs_0.25deg_regular_1979-2018_sub.nc"))

create.pet.makkink.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE)

############ PET Priestly-Taylor ############
pet.file <- "PET_priestley_testfile.nc"
out.file <- sprintf("%s%s", out.dir, pet.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"tn_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"tx_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"dp_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"rs_0.25deg_regular_1979-2018_sub.nc"))

create.pet.priestly.taylor.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE)

##### test output #######

## penman ss
in.nc_1= nc_open("PET_output/PET_penman_ss_testfile_master.nc")
lats <- ncvar_get( in.nc_1, "latitude")   # coordinate variable
nlat <- length(lats)
lons <- ncvar_get( in.nc_1, "longitude")   # coordinate variable
nlon <- length(lons)
tm <- ncvar_get( in.nc_1, "time")/24
nt <- length(tm)
indat1 <- as.matrix(ncvar_get( in.nc_1, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("PET_output/PET_penman_ss_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("penman_ss output matches master file")
}else{
  print("Warning: penman_ss output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)
## penman rs
in.nc_1= nc_open("PET_output/PET_penman_rs_testfile_master.nc")
indat1 <- as.matrix(ncvar_get( in.nc_1, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("PET_output/PET_penman_rs_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("penman_rs output matches master file")
}else{
  print("Warning: penman_rs output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)
## priestley
in.nc_1= nc_open("PET_output/PET_priestley_testfile_master.nc")
indat1 <- as.matrix(ncvar_get( in.nc_1, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("PET_output/PET_priestley_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("priestley output matches master file")
}else{
  print("Warning: priestley output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)
## makkink
in.nc_1= nc_open("PET_output/PET_makkink_testfile_master.nc")
indat1 <- as.matrix(ncvar_get( in.nc_1, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("PET_output/PET_makkink_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "pet", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("makkink output matches master file")
}else{
  print("Warning: makkink output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)