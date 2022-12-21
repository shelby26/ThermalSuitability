species <- "Acropora samoensis"
b80minas <- 15.7
b80maxas <- 26.4

species <- "Montipora caliculata"
b80minmc <- 16.1
b80maxmc <- 27.9

species <- "Galaxea fascicularis"
b80mingf <- 16.9
b80maxgf <- 27.1

library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plyr)
library(rworldmap)
require(sp)
require(raster)
require(rasterVis)

setwd("D:/FelixGIScomputerstuff/CMIP6-sst-climateprojections/cdotest/")


#1950-2014 11 models ensemble current - c
current<- nc_open("11models_historical_1850-2014.nc")

dname <- "tos"


lon <- ncvar_get(current,"lon")
nlonc <- dim(lon)

lat <- ncvar_get(current,"lat")
nlatc <- dim(lat)
timec <- ncvar_get(current,"time")

tunits <- ncatt_get(current,"time","units")

nt <- dim(timec)


tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
#chron(time,origin=c(tmonth, tday, tyear))


tos_array <- ncvar_get(current,dname)
dlname <- ncatt_get(current,dname,"long_name")
dunits <- ncatt_get(current,dname,"units")
fillvalue <- ncatt_get(current,dname,"_FillValue")

tos_array[tos_array==fillvalue$value] <- NA

##current decade slice
##using 30 year average, so 1984 to 2014
#january 1984 to 2014
m <- 1620:1980
tos_slice <- tos_array[,,m]
lonlat <- as.matrix(expand.grid(lon,lat))

tos_vec <- as.vector(tos_slice)

ntc <- 360
tos_mat <- matrix(tos_vec, nrow=nlonc*nlatc, ncol=ntc)

lonlat <- as.matrix(expand.grid(lon,lat))
tos_df01 <- data.frame(cbind(lonlat,tos_mat))


###species-specific b80 calculation
##proportion of time in b80 for monthly values for current temperatures

#Acropora samoensis
b80matc.as <- ifelse(tos_df01> b80minas & tos_df01< b80maxas, 1,0)
b80dfc.as <- data.frame(b80matc.as)
b80dfc.as$meanb80 <- apply(b80dfc.as[3:362],1,mean) 

#Montipora caliculata
b80matc.mc <- ifelse(tos_df01> b80minmc & tos_df01< b80maxmc, 1,0)
b80dfc.mc <- data.frame(b80matc.mc)
b80dfc.mc$meanb80 <- apply(b80dfc.mc[3:362],1,mean) 

#Galaxea fascicularis
b80matc.gf <- ifelse(tos_df01> b80mingf & tos_df01< b80maxgf, 1,0)
b80dfc.gf <- data.frame(b80matc.gf)
b80dfc.gf$meanb80 <- apply(b80dfc.gf[3:362],1,mean) 

##historical decade slice
##using 30 year average, so 1870 to 1900
#january 1870 to 1900
m <- 1:361
tos_slice <- tos_array[,,m]
lonlat <- as.matrix(expand.grid(lon,lat))

tos_vec <- as.vector(tos_slice)

nth <- 360
tos_mat <- matrix(tos_vec, nrow=nlonc*nlatc, ncol=nth)

lonlat <- as.matrix(expand.grid(lon,lat))
tos_df00 <- data.frame(cbind(lonlat,tos_mat))

###species-specific b80 calculation
##proportion of time in b80 for monthly values for historical

#Acropora samoensis
b80math.as <- ifelse(tos_df00> b80minas & tos_df00< b80maxas, 1,0)
b80dfh.as <- data.frame(b80math.as)
b80dfh.as$meanb80 <- apply(b80dfh.as[3:362],1,mean) 

#Montipora caliculata
b80math.mc <- ifelse(tos_df00> b80minmc & tos_df00< b80maxmc, 1,0)
b80dfh.mc <- data.frame(b80math.mc)
b80dfh.mc$meanb80 <- apply(b80dfh.mc[3:362],1,mean) 

#Galaxea fascicularis
b80math.gf <- ifelse(tos_df00> b80mingf & tos_df00< b80maxgf, 1,0)
b80dfh.gf <- data.frame(b80math.gf)
b80dfh.gf$meanb80 <- apply(b80dfh.gf[3:362],1,mean) 

##future 12 model ensemble decade slice
#2050-2100climate
future <- nc_open("12models_ssp585_2015-2100.nc")

dname <- "tos"

lonf <- ncvar_get(future,"lon")
nlonf <- dim(lon)

latf <- ncvar_get(future,"lat")
nlatf <- dim(lat)
timef <- ncvar_get(future,"time")

tunitsf <- ncatt_get(future,"time","units")

ntf <- dim(timef)


tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
#chron(time,origin=c(tmonth, tday, tyear))

tos_arrayf <- ncvar_get(future,dname)
dlnamef <- ncatt_get(future,dname,"long_name")
dunitsf <- ncatt_get(future,dname,"units")
fillvaluef <- ncatt_get(future,dname,"_FillValue")

tos_arrayf[tos_arrayf==fillvaluef$value] <- NA

#using 30 year average, so 2070 to 2100
n <- 672:1032
tos_slicef <- tos_arrayf[,,n]

lonlatf <- as.matrix(expand.grid(lonf,latf))

tos_vecf <- as.vector(tos_slicef)

ntf <- 360

tos_matf <- matrix(tos_vecf, nrow=nlonf*nlatf, ncol=ntf)

lonlatf <- as.matrix(expand.grid(lonf,latf))
tos_df02 <- data.frame(cbind(lonlatf,tos_matf))

###species-specific b80 calculation for future

#Acropora samoensis
b80matf.as <- ifelse(tos_df02> b80minas & tos_df02< b80maxas, 1,0)
b80dff.as <- data.frame(b80matf.as)
b80dff.as$meanb80 <- apply(b80dff.as[3:362],1,mean) 

#Montipora caliculata
b80matf.mc <- ifelse(tos_df02> b80minmc & tos_df02< b80maxmc, 1,0)
b80dff.mc <- data.frame(b80matf.mc)
b80dff.mc$meanb80 <- apply(b80dff.mc[3:362],1,mean) 

#Galaxea fascicularis
b80matf.gf <- ifelse(tos_df02> b80mingf & tos_df02< b80maxgf, 1,0)
b80dff.gf <- data.frame(b80matf.gf)
b80dff.gf$meanb80 <- apply(b80dff.gf[3:362],1,mean) 

##convert each into array, calculate delta

lon2 <- lonf
lat2 <- latf
time2 <- time
tunits2 <- tunits
nlon2 <- nlonf
nlat2 <- nlatf
nt2 <- nt

#Acropora samoensis

#current
mat2cb80.as <- as.matrix(b80dfc.as[3:(3+ntc-1)])
array2cb80.as  <- array(mat2cb80.as , dim=c(nlon2,nlat2,ntc))
array2c2b80.as <- array(b80dfc.as$meanb80, dim=c(nlon2,nlat2))

#historical
mat2hb80.as <- as.matrix(b80dfh.as[3:(3+nth-1)])
array2hb80.as  <- array(mat2hb80.as , dim=c(nlon2,nlat2,nth))
array2h2b80.as <- array(b80dfh.as$meanb80, dim=c(nlon2,nlat2))

#future
mat2fb80.as <- as.matrix(b80dff.as[3:(3+ntf-1)])
array2fb80.as  <- array(mat2fb80.as , dim=c(nlon2,nlat2,ntf))
array2f2b80.as <- array(b80dff.as$meanb80, dim=c(nlon2,nlat2))

#Montipora caliculata

#current
mat2cb80.mc <- as.matrix(b80dfc.mc[3:(3+ntc-1)])
array2cb80.mc  <- array(mat2cb80.mc , dim=c(nlon2,nlat2,ntc))
array2c2b80.mc <- array(b80dfc.mc$meanb80, dim=c(nlon2,nlat2))

#historical
mat2hb80.mc <- as.matrix(b80dfh.mc[3:(3+nth-1)])
array2hb80.mc  <- array(mat2hb80.mc , dim=c(nlon2,nlat2,ntc))
array2h2b80.mc <- array(b80dfh.mc$meanb80, dim=c(nlon2,nlat2))

#future
mat2fb80.mc <- as.matrix(b80dff.mc[3:(3+ntf-1)])
array2fb80.mc  <- array(mat2fb80.mc , dim=c(nlon2,nlat2,ntf))
array2f2b80.mc <- array(b80dff.mc$meanb80, dim=c(nlon2,nlat2))

#Galaxea fascicularis

#current
mat2cb80.gf <- as.matrix(b80dfc.gf[3:(3+ntc-1)])
array2cb80.gf  <- array(mat2cb80.gf , dim=c(nlon2,nlat2,ntc))
array2c2b80.gf <- array(b80dfc.gf$meanb80, dim=c(nlon2,nlat2))

#historical
mat2hb80.gf <- as.matrix(b80dfh.gf[3:(3+nth-1)])
array2hb80.gf  <- array(mat2hb80.gf , dim=c(nlon2,nlat2,nth))
array2h2b80.gf <- array(b80dfh.gf$meanb80, dim=c(nlon2,nlat2))

#future
mat2fb80.gf <- as.matrix(b80dff.gf[3:(3+ntf-1)])
array2fb80.gf  <- array(mat2fb80.gf , dim=c(nlon2,nlat2,ntf))
array2f2b80.gf <- array(b80dff.gf$meanb80, dim=c(nlon2,nlat2))

#change over time
#Acropora samoensis
#current to future
delta_array.as <- array2f2b80.as - array2c2b80.as
#historical to current
delta_array1.as <- array2c2b80.as - array2h2b80.as

#Montipora caliculata
#current to future
delta_array.mc <- array2f2b80.mc - array2c2b80.mc
#historical to current
delta_array1.mc <- array2c2b80.mc - array2h2b80.mc

#Galaxea fascicularis
#current to future
delta_array.gf <- array2f2b80.gf - array2c2b80.gf
#historical to current
delta_array1.gf <- array2c2b80.gf - array2h2b80.gf


#plotting change

grid <- expand.grid(lon=lon, lat=lat)

data("coastsCoarse")

cutpts<- seq(-1,1,0.05)
mapTheme <- rasterTheme(region = brewer.pal(10, "RdBu"))

#Acropora samoensis
#historical to current
as.h.to.c <- levelplot(delta_array1.as ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('AS')))) + layer(sp.polygons(coastsCoarse))
#current to future
as.c.to.f <- levelplot(delta_array.as ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('AS')))) + layer(sp.polygons(coastsCoarse))

#Montipora caliculata
#historical to current
mc.h.to.c <- levelplot(delta_array1.mc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('MC')))) + layer(sp.polygons(coastsCoarse))
#current to future
mc.c.to.f <- levelplot(delta_array.mc ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('MC')))) + layer(sp.polygons(coastsCoarse))

#Galaxea fascicularis
#historical to current
gf.h.to.c <- levelplot(delta_array1.gf ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('GF')))) + layer(sp.polygons(coastsCoarse))
#current to future
gf.c.to.f <- levelplot(delta_array.gf ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('GF')))) + layer(sp.polygons(coastsCoarse))


tiff("as.h.to.c.tiff", units="in", width=7, height=5, res=600)
as.h.to.c
dev.off()

tiff("as.c.to.f.tiff", units="in", width=7, height=5, res=600)
as.c.to.f
dev.off()

tiff("mc.h.to.c.tiff", units="in", width=7, height=5, res=600)
mc.h.to.c
dev.off()

tiff("mc.c.to.f.tiff", units="in", width=7, height=5, res=600)
mc.c.to.f
dev.off()

tiff("gf.h.to.c.tiff", units="in", width=7, height=5, res=600)
gf.h.to.c
dev.off()

tiff("gf.c.to.f.tiff", units="in", width=7, height=5, res=600)
gf.c.to.f
dev.off()


library(ggpubr)
ggarrange()
