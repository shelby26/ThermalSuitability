
##data list
setwd("F:\\FME2020_daily\\FME_data_copy\\Tanaka2019_SICAT02_MIROC_MRI\\tanakay01.yes.jamstec.go.jp\\Data\\SICAT02\\MRI-CGCM3\\result_all")

present <- dir(,"MRI-CGCM3_rcp85_to-sfc_.*_mean.tif")[1:120]
near2050 <- dir(,"MRI-CGCM3_rcp85_to-sfc_.*_mean.tif")[121:252]
far2100 <- dir(,"MRI-CGCM3_rcp85_to-sfc_.*_mean.tif")[253:432]

require(raster)

##data extraction
#filenames <- present
main <- function(filenames){
	i <- 1
	 r <- raster(filenames[i])
	 df1 <- as.data.frame(r, xy=TRUE)
	for(i in 2:length(filenames)){
	 r <- raster(filenames[i])
	 df1 <- cbind(df1,as.data.frame(r))
	}
	return(df1)
}

tos_df001.P <- main(present)
tos_df001.N <- main(near2050)
tos_df001.F <- main(far2100)



############
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


###species-specific b80 calculation
##proportion of time in b80 for monthly values for ""Present""

#Oulastrea crispata
PR2matP.oc <- ifelse(tos_df001.P> PR2minoc, 1,0)
PR2dfP.oc <- data.frame(PR2matP.oc)
PR2dfP.oc$meanPR2 <- apply(PR2dfP.oc[-(1:2)],1,mean)

#Porites lobata
PR2matP.pl <- ifelse(tos_df001.P> PR2minpl, 1,0)
PR2dfP.pl <- data.frame(PR2matP.pl)
PR2dfP.pl$meanPR2 <- apply(PR2dfP.pl[-(1:2)],1,mean)


###species-specific b80 calculation for ""Nearfuture""

#Oulastrea crispata
PR2matN.oc <- ifelse(tos_df001.N> PR2minoc, 1,0)
PR2dfN.oc <- data.frame(PR2matN.oc)
PR2dfN.oc$meanPR2 <- apply(PR2dfN.oc[-(1:2)],1,mean)

#Porites lobata
PR2matN.pl <- ifelse(tos_df001.N> PR2minpl, 1,0)
PR2dfN.pl <- data.frame(PR2matN.pl)
PR2dfN.pl$meanPR2 <- apply(PR2dfN.pl[-(1:2)],1,mean)

###species-specific b80 calculation for ""Farfuture""

#Oulastrea crispata
PR2matF.oc <- ifelse(tos_df001.F> PR2minoc, 1,0)
PR2dfF.oc <- data.frame(PR2matF.oc)
PR2dfF.oc$meanPR2 <- apply(PR2dfF.oc[-(1:2)],1,mean)

#Porites lobata
PR2matF.pl <- ifelse(tos_df001.F> PR2minpl, 1,0)
PR2dfF.pl <- data.frame(PR2matF.pl)
PR2dfF.pl$meanPR2 <- apply(PR2dfF.pl[-(1:2)],1,mean)



##
##convert each into array, calculate delta
##convert each into array, calculate delta

lon2 <- tos_df001.P[,1]
lat2 <- tos_df001.P[,2]
#time2 <- time
#tunits2 <- tunits
nlon2 <- length(lon2)
nlat2 <- length(lat2)
#nt2 <- nt

ntcP <- length(present)

save(list=ls(), file="all.Rdata") #

####<<<latter comment outed ways were not work because of the memory size.>>>
#did you just wanted to see annual change of the image? if so I recommend to use rasterFromXYZ().
#
###Oulastrea crispata
##
###Present
##mat2cPR2.oc <- as.matrix(PR2dfP.oc[,-(1:2)][,1:ntcP])        #I think we need connma? not PR2dfP.oc[-(1:2)] but PR2dfP.oc[,-(1:2)]
##array2cPR2.oc  <- array(mat2cPR2.oc , dim=c(nlon2,nlat2,ntcP))
##array2c2PR2.oc <- array(PR2dfP.oc$meanPR2, dim=c(nlon2,nlat2))
##
###Nearfuture
##mat2hPR2.oc <- as.matrix(PR2dfN.oc[-(1:2)])
##array2hPR2.oc  <- array(mat2hPR2.oc , dim=c(nlon2,nlat2,ntc))
##array2h2PR2.oc <- array(PR2dfN.oc$meanPR2, dim=c(nlon2,nlat2))
##
###Farfuture
##mat2fPR2.oc <- as.matrix(PR2dfF.oc[-(1:2)])
##array2fPR2.oc  <- array(mat2fPR2.oc , dim=c(nlon2,nlat2,ntc))
##array2f2PR2.oc <- array(PR2dfF.oc$meanPR2, dim=c(nlon2,nlat2))
##

##<<<so other ways only for the mean values are culculated>>>
#
xyz.list <- list()
xyz.list[[1]] <- data.frame(lon2,lat2,PR2dfP.oc$meanPR2)
xyz.list[[2]] <- data.frame(lon2,lat2,PR2dfN.oc$meanPR2)
xyz.list[[3]] <- data.frame(lon2,lat2,PR2dfF.oc$meanPR2)

xyz.list[[4]] <- data.frame(lon2,lat2,PR2dfP.pl$meanPR2)
xyz.list[[5]] <- data.frame(lon2,lat2,PR2dfN.pl$meanPR2)
xyz.list[[6]] <- data.frame(lon2,lat2,PR2dfF.pl$meanPR2)

#<<<change over time only for the mean calues is this diffrent from orijinal?>>>
#Oulastrea crispata
#present to nearfuture
delta_array001N_P.oc <- PR2dfN.oc$meanPR2 - PR2dfP.oc$meanPR2
#present to farfuture
delta_array001F_P.oc <- PR2dfF.oc$meanPR2 - PR2dfP.oc$meanPR2

#Porites lobata
#present to nearfuture
delta_array001N_P.pl <- PR2dfN.pl$meanPR2 - PR2dfP.pl$meanPR2
#present to farfuture
delta_array001F_P.pl <- PR2dfF.pl$meanPR2 - PR2dfP.pl$meanPR2


xyz.list[[ 7]] <- data.frame(lon2,lat2,delta_array001N_P.oc)
xyz.list[[ 8]] <- data.frame(lon2,lat2,delta_array001F_P.oc)
xyz.list[[ 9]] <- data.frame(lon2,lat2,delta_array001N_P.pl)
xyz.list[[10]] <- data.frame(lon2,lat2,delta_array001F_P.pl)

save(list=xyz.list, file="xyz.list.Rdata") #<<<light data for plotting>>>

###########################################
##for drowing data

#<<<plot images>>>
load("xyz.list.Rdata") 

ras.test1 <- rasterFromXYZ(xyz.list[[1]], res=c(NA,NA), crs="", digits=5)
image(ras.test1)

library(rgdal)
install.packages("MetBrewer")
install.packages("devtools") 
devtools::install_github("BlakeRMills/MetBrewer") 
library(MetBrewer)
library(rasterVis)
library(RColorBrewer)


cutpts<- seq(-1,1,0.05)
cutpts1<- seq(0,1,0.01)
nice <- rasterTheme(region = met.brewer("Hokusai1",n=100))   #I think Hokusai is not work and it mean hokusai1 or hokusai2 ?
mapTheme <- rasterTheme(region = brewer.pal(10, "RdBu"))
ylgnbu <- rasterTheme(region = brewer.pal(7, "YlGnBu"))


lon2 <- xyz.list[[1]][,1]
lat2 <- xyz.list[[1]][,2]
val <- xyz.list[[1]][,3]

#test plot of each result
levelplot(val ~ lon2 * lat2,  colorkey=TRUE,    ,    , at=cutpts1, pretty=T, 
          par.settings = ylgnbu, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC'))))


#Present to near future
levelplot(xyz.list[[7]][,3]~ lon2 * lat2,  colorkey=TRUE,      ,      , at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC'))))
#present to far future
levelplot(xyz.list[[8]][,3]~ lon2 * lat2,  colorkey=TRUE,      ,      , at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('OC'))))
#present to near future
levelplot(xyz.list[[9]][,3]~ lon2 * lat2,  colorkey=TRUE,      ,      , at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('PL'))))
#present to far future
levelplot(xyz.list[[10]][,3]~ lon2 * lat2,  colorkey=TRUE,      ,      , at=cutpts, pretty=T, 
                       par.settings = mapTheme, xlab="Longitude",  ylab="Latitude", main= substitute(paste(italic('PL'))))

#
#
#######
