species <- "Porites lobata"
PR2minpl <- 9.3
PR2maxpl <- 32

species <- "Oulastrea crispata"
PR2minoc <- 18.1
PR2maxoc <- 32 

species <- "Acropora samoensis"
PR2minas <- 8.4
PR2maxas <- 33.7

species <- "Galaxea fascicularis"
PR2mingf <- 9.2
PR2maxgf <- 31.1

species <- "Montipora caliculata"
PR2minmc <- 11.3
PR2maxmc <- 32.6 



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

b <- brick("11models_historical_1850-2014.nc", varname="tos")

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

#Oulastrea crispata
PR2matc.oc <- ifelse(tos_df01> PR2minoc, 1,0)
PR2dfc.oc <- data.frame(PR2matc.oc)
PR2dfc.oc$meanPR2 <- apply(PR2dfc.oc[3:362],1,mean)

#Porites lobata
PR2matc.pl <- ifelse(tos_df01> PR2minpl, 1,0)
PR2dfc.pl <- data.frame(PR2matc.pl)
PR2dfc.pl$meanPR2 <- apply(PR2dfc.pl[3:362],1,mean)

#Acropora samoensis
PR2matc.as <- ifelse(tos_df01> PR2minas, 1,0)
PR2dfc.as <- data.frame(PR2matc.as)
PR2dfc.as$meanPR2 <- apply(PR2dfc.as[3:362],1,mean)

#Galaxea fascicularis
PR2matc.gf <- ifelse(tos_df01> PR2mingf, 1,0)
PR2dfc.gf <- data.frame(PR2matc.gf)
PR2dfc.gf$meanPR2 <- apply(PR2dfc.gf[3:362],1,mean)

#Montipora caliculata
PR2matc.mc <- ifelse(tos_df01> PR2minmc, 1,0)
PR2dfc.mc <- data.frame(PR2matc.mc)
PR2dfc.mc$meanPR2 <- apply(PR2dfc.mc[3:362],1,mean)


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

#Oulastrea crispata
PR2math.oc <- ifelse(tos_df00> PR2minoc, 1,0)
PR2dfh.oc <- data.frame(PR2math.oc)
PR2dfh.oc$meanPR2 <- apply(PR2dfh.oc[3:362],1,mean)

#Porites lobata
PR2math.pl <- ifelse(tos_df00> PR2minpl, 1,0)
PR2dfh.pl <- data.frame(PR2math.pl)
PR2dfh.pl$meanPR2 <- apply(PR2dfh.pl[3:362],1,mean)

#Acropora samoensis
PR2math.as <- ifelse(tos_df00> PR2minas, 1,0)
PR2dfh.as <- data.frame(PR2math.as)
PR2dfh.as$meanPR2 <- apply(PR2dfh.as[3:362],1,mean)

#Galaxea fascicularis
PR2math.gf <- ifelse(tos_df00> PR2mingf, 1,0)
PR2dfh.gf <- data.frame(PR2math.gf)
PR2dfh.gf$meanPR2 <- apply(PR2dfh.gf[3:362],1,mean)

#Montipora caliculata
PR2math.mc <- ifelse(tos_df00> PR2minmc, 1,0)
PR2dfh.mc <- data.frame(PR2math.mc)
PR2dfh.mc$meanPR2 <- apply(PR2dfh.mc[3:362],1,mean)
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

#Oulastrea crispata
PR2matf.oc <- ifelse(tos_df02> PR2minoc, 1,0)
PR2dff.oc <- data.frame(PR2matf.oc)
PR2dff.oc$meanPR2 <- apply(PR2dff.oc[3:362],1,mean)

#Porites lobata
PR2matf.pl <- ifelse(tos_df02> PR2minpl, 1,0)
PR2dff.pl <- data.frame(PR2matf.pl)
PR2dff.pl$meanPR2 <- apply(PR2dff.pl[3:362],1,mean)

#Acropora samoensis
PR2matf.as <- ifelse(tos_df02> PR2minas, 1,0)
PR2dff.as <- data.frame(PR2matf.as)
PR2dff.as$meanPR2 <- apply(PR2dff.as[3:362],1,mean)

#Galaxea fascicularis
PR2matf.gf <- ifelse(tos_df02> PR2mingf, 1,0)
PR2dff.gf <- data.frame(PR2matf.gf)
PR2dff.gf$meanPR2 <- apply(PR2dff.gf[3:362],1,mean)

#Montipora caliculata
PR2matf.mc <- ifelse(tos_df02> PR2minmc, 1,0)
PR2dff.mc <- data.frame(PR2matf.mc)
PR2dff.mc$meanPR2 <- apply(PR2dff.mc[3:362],1,mean)


##convert each into array, calculate delta

lon2 <- lonf
lat2 <- latf
time2 <- time
tunits2 <- tunits
nlon2 <- nlonf
nlat2 <- nlatf
nt2 <- nt

#Oulastrea crispata

#current
mat2cPR2.oc <- as.matrix(PR2dfc.oc[3:(3+ntc-1)])
array2cPR2.oc  <- array(mat2cPR2.oc , dim=c(nlon2,nlat2,ntc))
array2c2PR2.oc <- array(PR2dfc.oc$meanPR2, dim=c(nlon2,nlat2))

#historical
mat2hPR2.oc <- as.matrix(PR2dfh.oc[3:(3+ntc-1)])
array2hPR2.oc  <- array(mat2hPR2.oc , dim=c(nlon2,nlat2,ntc))
array2h2PR2.oc <- array(PR2dfh.oc$meanPR2, dim=c(nlon2,nlat2))

#future
mat2fPR2.oc <- as.matrix(PR2dff.oc[3:(3+ntc-1)])
array2fPR2.oc  <- array(mat2fPR2.oc , dim=c(nlon2,nlat2,ntc))
array2f2PR2.oc <- array(PR2dff.oc$meanPR2, dim=c(nlon2,nlat2))

#Porites lobata

#current
mat2cPR2.pl <- as.matrix(PR2dfc.pl[3:(3+ntc-1)])
array2cPR2.pl  <- array(mat2cPR2.pl , dim=c(nlon2,nlat2,ntc))
array2c2PR2.pl <- array(PR2dfc.pl$meanPR2, dim=c(nlon2,nlat2))

#historical
mat2hPR2.pl <- as.matrix(PR2dfh.pl[3:(3+ntc-1)])
array2hPR2.pl  <- array(mat2hPR2.pl , dim=c(nlon2,nlat2,ntc))
array2h2PR2.pl <- array(PR2dfh.pl$meanPR2, dim=c(nlon2,nlat2))

#future
mat2fPR2.pl <- as.matrix(PR2dff.pl[3:(3+ntc-1)])
array2fPR2.pl  <- array(mat2fPR2.pl , dim=c(nlon2,nlat2,ntc))
array2f2PR2.pl <- array(PR2dff.pl$meanPR2, dim=c(nlon2,nlat2))

#Acropora samoensis

#current
mat2cPR2.as <- as.matrix(PR2dfc.as[3:(3+ntc-1)])
array2cPR2.as  <- array(mat2cPR2.as , dim=c(nlon2,nlat2,ntc))
array2c2PR2.as <- array(PR2dfc.as$meanPR2, dim=c(nlon2,nlat2))

#historical
mat2hPR2.as <- as.matrix(PR2dfh.as[3:(3+ntc-1)])
array2hPR2.as  <- array(mat2hPR2.as , dim=c(nlon2,nlat2,ntc))
array2h2PR2.as <- array(PR2dfh.as$meanPR2, dim=c(nlon2,nlat2))

#future
mat2fPR2.as <- as.matrix(PR2dff.as[3:(3+ntc-1)])
array2fPR2.as  <- array(mat2fPR2.as , dim=c(nlon2,nlat2,ntc))
array2f2PR2.as <- array(PR2dff.as$meanPR2, dim=c(nlon2,nlat2))

#Galaxea fascicularis

#current
mat2cPR2.gf <- as.matrix(PR2dfc.gf[3:(3+ntc-1)])
array2cPR2.gf  <- array(mat2cPR2.gf , dim=c(nlon2,nlat2,ntc))
array2c2PR2.gf <- array(PR2dfc.gf$meanPR2, dim=c(nlon2,nlat2))

#historical
mat2hPR2.gf <- as.matrix(PR2dfh.gf[3:(3+ntc-1)])
array2hPR2.gf  <- array(mat2hPR2.gf , dim=c(nlon2,nlat2,ntc))
array2h2PR2.gf <- array(PR2dfh.gf$meanPR2, dim=c(nlon2,nlat2))

#future
mat2fPR2.gf <- as.matrix(PR2dff.gf[3:(3+ntc-1)])
array2fPR2.gf  <- array(mat2fPR2.gf , dim=c(nlon2,nlat2,ntc))
array2f2PR2.gf <- array(PR2dff.gf$meanPR2, dim=c(nlon2,nlat2))


#Montipora caliculata

#current
mat2cPR2.mc <- as.matrix(PR2dfc.mc[3:(3+ntc-1)])
array2cPR2.mc  <- array(mat2cPR2.mc , dim=c(nlon2,nlat2,ntc))
array2c2PR2.mc <- array(PR2dfc.mc$meanPR2, dim=c(nlon2,nlat2))

#historical
mat2hPR2.mc <- as.matrix(PR2dfh.mc[3:(3+ntc-1)])
array2hPR2.mc  <- array(mat2hPR2.mc , dim=c(nlon2,nlat2,ntc))
array2h2PR2.mc <- array(PR2dfh.mc$meanPR2, dim=c(nlon2,nlat2))

#future
mat2fPR2.mc <- as.matrix(PR2dff.mc[3:(3+ntc-1)])
array2fPR2.mc  <- array(mat2fPR2.mc , dim=c(nlon2,nlat2,ntc))
array2f2PR2.mc <- array(PR2dff.mc$meanPR2, dim=c(nlon2,nlat2))

#allocate extra memory
memory.limit()
memory.limit(size=18000)

#change over time
#Oulastrea crispata
#current to future
delta_array.oc <- array2f2PR2.oc - array2c2PR2.oc
#historical to current
delta_array1.oc <- array2c2PR2.oc - array2h2PR2.oc

#Porites lobata
#current to future
delta_array.pl <- array2f2PR2.pl - array2c2PR2.pl
#historical to current
delta_array1.pl <- array2c2PR2.pl - array2h2PR2.pl

#Acropora samoensis
#current to future
delta_array.as <- array2f2PR2.as - array2c2PR2.as
#historical to current
delta_array1.as <- array2c2PR2.as - array2h2PR2.as

#Galaxea fascicularis
#current to future
delta_array.gf <- array2f2PR2.gf - array2c2PR2.gf
#historical to current
delta_array1.gf <- array2c2PR2.gf - array2h2PR2.gf

#Montipora caliculata
#current to future
delta_array.mc <- array2f2PR2.mc - array2c2PR2.mc
#historical to current
delta_array1.mc <- array2c2PR2.mc - array2h2PR2.mc


grid <- expand.grid(lon=lon, lat=lat)

data("coastsCoarse")
library(rgdal)
library(MetBrewer)
install.packages("MetBrewer")

install.packages("devtools") 
devtools::install_github("BlakeRMills/MetBrewer")

coast<- readOGR(file.choose()) #Coastline shapefile GSHHS_l_L1.shp for "low" res

cutpts<- seq(-1,1,0.05)
cutpts1<- seq(0,1,0.01)
nice <- rasterTheme(region = met.brewer("Cassatt1",n=100))
mapTheme <- rasterTheme(region = brewer.pal(10, "RdBu"))
ylgnbu <- rasterTheme(region = brewer.pal(9, "PuBuGn"))



levelplot(array2h2PR2.oc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(100,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coast, size = 1))

levelplot(array2c2PR2.oc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coastsCoarse))

levelplot(array2f2PR2.oc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coastsCoarse))


levelplot(array2h2PR2.pl ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coastsCoarse))

levelplot(array2c2PR2.pl ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coastsCoarse))

levelplot(array2f2PR2.pl ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
          par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coastsCoarse))

#Oulastrea crispata
#historical to current
oc.h.to.c <- levelplot(delta_array1.oc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coast))
#current to future
oc.c.to.f <- levelplot(delta_array.oc ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC')))) + layer(sp.polygons(coast))

#Porites lobata
#historical to current
pl.h.to.c <- levelplot(delta_array1.pl ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('PL')))) + layer(sp.polygons(coast))
#current to future
pl.c.to.f <- levelplot(delta_array.pl ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('PL')))) + layer(sp.polygons(coast))

#Acropora samoensis
#historical to current
as.h.to.c <- levelplot(delta_array1.as ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('AS')))) + layer(sp.polygons(coast))
#current to future
as.c.to.f <- levelplot(delta_array.as ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('AS')))) + layer(sp.polygons(coast))

#Galaxea fascicularis
#historical to current
gf.h.to.c <- levelplot(delta_array1.gf ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('GF')))) + layer(sp.polygons(coast))
#current to future
gf.c.to.f <- levelplot(delta_array.gf ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('GF')))) + layer(sp.polygons(coast))

#Montipora caliculata
mc.h.to.c <- levelplot(delta_array1.mc ~ lon * lat,  colorkey=TRUE,data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings =ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('MC')))) + layer(sp.polygons(coast))
#current to future
mc.c.to.f <- levelplot(delta_array.mc ~ lon * lat, colorkey=TRUE, data=grid, xlim=c(90,170), ylim=c(10,60), at=cutpts1, pretty=T, 
                       par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('MC')))) + layer(sp.polygons(coast))


tiff("pr2.oc.h.to.c.tiff", units="in", width=7, height=5, res=600)
oc.h.to.c
dev.off()

tiff("pr2.oc.c.to.f.tiff", units="in", width=7, height=5, res=600)
oc.c.to.f
dev.off()


tiff("pr2.pl.h.to.c.tiff", units="in", width=7, height=5, res=600)
pl.h.to.c
dev.off()

tiff("pr2.pl.c.to.f.tiff", units="in", width=7, height=5, res=600)
pl.c.to.f
dev.off()


tiff("pr2.as.h.to.c.tiff", units="in", width=7, height=5, res=600)
as.h.to.c
dev.off()

tiff("pr2.as.c.to.f.tiff", units="in", width=7, height=5, res=600)
as.c.to.f
dev.off()


tiff("pr2.gf.h.to.c.tiff", units="in", width=7, height=5, res=600)
gf.h.to.c
dev.off()

tiff("pr2.gf.c.to.f.tiff", units="in", width=7, height=5, res=600)
gf.c.to.f
dev.off()

tiff("pr2.mc.h.to.c.tiff", units="in", width=7, height=5, res=600)
mc.h.to.c
dev.off()

tiff("pr2.mc.c.to.f.tiff", units="in", width=7, height=5, res=600)
mc.c.to.f
dev.off()