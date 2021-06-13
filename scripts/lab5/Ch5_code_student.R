########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 5: spatial dependence
########################################################
########################################################

#load packages
library(pgirmess)         #for simple correlograms; version 1.6.9 used
library(ncf)              #for spline correlograms; version 1.2-5 used
library(spdep)            #for general correlograms; version 0.7-8 used
library(geoR)             #for variograms/kriging; version 1.7-5.2.1 used
library(gstat)            #for variograms/kriging; version 1.1-6 used
library(RandomFields)     #for simulating spatial random fields/spatial dependence; version 3.1.50 used
library(raster)           #for raster covariate data; version 2.6-7 used
library(waveslim)         #for multi-scale dependence; version 1.7.5 used
library(fields)           #for multi-scale dependence; version 9.6 used
library(reshape2)         #for re-formatting data; version 1.4.3 used
library(vegan)            #for multi-scale dependence; version 2.5-2 used
library(adespatial)       #for multi-scale dependence; version 0.3-0 used

#set working directory where data were downloaded
setwd("~/Desktop/GitHub/ECOL_620/data/data_for_lab5")

###############################################
#5.3.3 Correlograms
###############################################

#load data
matrix <- read.csv('cactus_matrix.csv', header=T)

#inspect
head(matrix)

#plot
plot(matrix[,"y"] ~ matrix[,"x"],
     pch=21, cex=1.2,
     bg=gray.colors(12)[cut(matrix[,"Height"], breaks = 12)])

#plot distribution of vegetation height, not normally dist 
hist(matrix[,"Height"], xlab="Vegetation height (cm)",
     main="histogram of vegetation height")

#calculate distance matrix
coords <- cbind(matrix$x, matrix$y)
colnames(coords) <- c("x", "y")
distmat <- as.matrix(dist(coords))

#inspect
dim(distmat)
max(distmat)

#Maximum distance to consider in correlogram/variogram (~1/2 to 2/3 total dist) don't go beyond 47-50 m
maxdist <- 2/3*max(distmat)
maxdist

#--------------------------------------#
#Correlogram with pgirmess
#--------------------------------------#
# using height, moran i, nbclasses is related to pairwise comparisons 
correlog.pgirmess <- pgirmess::correlog(coords, matrix$Height, method="Moran",
                            nbclass=14, alternative = "two.sided")

#Moran and P values for each distance class, coef is moran i, pval- is moran i sig diff than 0, n is number pointwise comparisons in that distance band
round(correlog.pgirmess,2)

#plot: high correlation at 0.19, starts to break down and after 50 m p bad 
plot(correlog.pgirmess[,1], correlog.pgirmess[,2],
     xlab="Distance (m)", ylab="Moran's I", col = ifelse(correlog.pgirmess[,3] < 0.05,'red','black'), pch=19)
abline(h=0)

#--------------------------------------#
#Correlograms with ncf 
#--------------------------------------#

#Correlogram with non-parameteric test of significance, spline and conf int
correlog.ncf <- ncf::correlog(x = matrix$x, y = matrix$y, z = matrix$Height,
                         increment=5, resamp=99)

#plot
plot(correlog.ncf)
abline(h=0)

#Spline correlogram with 95% pointwise bootstrap confidence intervals
spline.corr <- spline.correlog(x = matrix$x, y = matrix$y, z = matrix$Height,
                               xmax = maxdist, resamp=99, type="boot")

#plot
plot(spline.corr)

#--------------------------------------#
#Moran's I test with spdep
#--------------------------------------#
# pull out the relevant stuff 
#make a neighborhood list
neigh <- dnearneigh(x=coords, d1=0, d2=3, longlat=F)#d1 is minimum distance, d2 is max distance

#plot the neighorhood
plot(neigh,coordinates(coords))

#create weights for the neighbors
wts <- nb2listw(neighbours=neigh, style='W', zero.policy=T)#W = row-standardized weights

#Moran's I test with normal approximation versus Monte Carlo permutation test
mor.mc <- moran.mc(x=matrix$Height, listw=wts, nsim=999, zero.policy=T) #Monte Carlo
mor.norm <- moran.test(x=matrix$Height, listw=wts, randomisation=F, zero.policy=T)#normal approximation

#inspect
mor.mc
mor.norm

#--------------------------------------#
#Correlograms with spdep
#--------------------------------------#

#repeat above with a for loop using consecutive lag distances

#Create data frame for storing output
correlog.sp <- data.frame(dist=seq(5, maxdist, by=5),
                        Morans.i=NA, Null.lcl=NA, Null.ucl=NA, Pvalue=NA)
#inspect
head(correlog.sp)

#then do a for loop to calculate Moran's I for lag distances
for (i in 1:nrow(correlog.sp)){

  d.start <- correlog.sp[i,"dist"]-5
  d.end <- correlog.sp[i,"dist"]

  neigh <- dnearneigh(x=coords, d1=d.start, d.end, longlat=F)
  wts <- nb2listw(neighbours=neigh, style='W', zero.policy=T)
  mor.i <- moran.mc(x=matrix$Height, listw=wts, nsim=99, alternative="greater", zero.policy=T)

  #summarize results from spdep
  correlog.sp[i, "dist"] <- (d.end+d.start)/2                                    #mean dist
  correlog.sp[i, "Morans.i"] <- mor.i$statistic 								                 #observed Moran's I
  correlog.sp[i, "Null.lcl"] <- quantile(mor.i$res, probs = 0.025,na.rm = TRUE)  #lower null envelope
  correlog.sp[i, "Null.ucl"] <- quantile(mor.i$res, probs = 0.975,na.rm = TRUE)  #upper null envelope
  correlog.sp[i, "Pvalue"] <- mor.i$p.value									                     #p-value for Moran's I at that distance category
}

#plot
plot(y=correlog.sp$Morans.i, x=correlog.sp$dist,
     xlab="Lag Distance(m)", ylab="Moran's I", ylim=c(-0.3,0.3))         #ylim provides limit on y-axis between -1 and 1
abline(h=0)                                                              #0 reference
lines(correlog.sp$dist, correlog.sp$Null.lcl,col = "red")	               #add the null lcl to the plot
lines(correlog.sp$dist, correlog.sp$Null.ucl,col = "red")	               #add the null ucl to the plot


#############################################
#5.3.4 Variograms
#############################################

#----------------------------------#
#in geoR
#----------------------------------#

#create geoR object
geoR.veg <- as.geodata(matrix)

#plot
plot(geoR.veg)

#Empirical semivariogram
emp <- variog(geoR.veg, max.dist=maxdist)

#plot, y is semivariance, plateau 
plot(emp)

#standardize breaks
emp <- variog(geoR.veg, max.dist=maxdist, breaks=c(seq(0,maxdist,by=3)))

#plot variogram
plot(emp)

#----------------------------------#
#in gstat
#----------------------------------#

gstat.veg <- matrix
coordinates(gstat.veg) = ~x + y

#Empirical semivariogram
emp.gstat <- variogram(Height ~ 1, cutoff=maxdist, width=3, gstat.veg)

#plot variogram
plot(emp.gstat)

#----------------------------------#
#anisotropy check
#----------------------------------#

#Directional variogram in geoR
emp4 <- variog4(geoR.veg, max.dist=maxdist)

#plot directional variogram
plot(emp4)

#Directional variogram in gstat
emp4.gstat <- variogram(Height ~ 1, cutoff=maxdist, alpha=c(0,45,90,135), gstat.veg)

#plot directional variogram
plot(emp4.gstat)

#-------------------------------------------------------#
#Model-based variograms with likelihood fitting in geoR
#-------------------------------------------------------#

#Exponential variogram
mlexp <- likfit(geoR.veg, cov.model="exp", ini=c(700,10))

#Spherical variogram
mlsph <- likfit(geoR.veg, cov.model="sph", ini=c(700,10))

#inspect
summary(mlexp)
summary(mlsph)
AIC(mlexp,mlsph) #model selection tool to look for models w lower aic values 

#plot
plot(emp)
lines(mlexp, col="blue") #note this appears as a poorer fit, as compared to "mlsph". But extracting values and plotting in Fig. 5.7 shows better fit
lines(mlsph, col="red")

#Monte Carlo envelopes
emp.env <- variog.mc.env(geoR.veg, obj.var=emp)

#plot
plot(emp, envelope=emp.env)
lines(mlsph, col="red")
lines(mlexp, col="blue")

#--------------------------------------------------#
#Model-based variograms with least squares in gstat
#--------------------------------------------------#

#inspect models that gstat can fit
vgm()#fits more models than geoR
show.vgms()#plots examples of the various variograms

#Spherical variogram
sph.gstat <- fit.variogram(emp.gstat, vgm("Sph")) #in vgm(psill, model, range, nugget)

#Exponential variogram
exp.gstat <- fit.variogram(emp.gstat, vgm("Exp"))

#inspect
exp.gstat
sph.gstat

#plot
plot(emp.gstat, exp.gstat)
plot(emp.gstat, sph.gstat)
  
#################################################
#5.3.5 Kriging and interpolation
#################################################

#Create grid with intervals of 0.5 unit (.5-m)
new.grid.5m <- expand.grid(x=seq(0,max(matrix$x),.5), y=seq(0,max(matrix$y),.5))#need labels for coords

#Maximum likelihood estimates from the variogram modeling
mlexp$nugget # nugget
mlexp$cov.pars[1] #partial sill
mlexp$cov.pars[2] #range

#Ordinary kriging
krig.geoR.exp <- krige.conv(geoR.veg,locations=new.grid.5m,
                          krige=krige.control(cov.pars=c(mlexp$cov.pars[1], mlexp$cov.pars[2]), nugget = mlexp$nugget,
                                              cov.model="exp", type.krige="OK"))

#inspect
hist(krig.geoR.exp$predict)
hist(matrix$Height)

#plot
image(krig.geoR.exp, main="krigged estimates") #see where points were collected 
image(krig.geoR.exp, val=sqrt(krig.geoR.exp$krige.var), main="kriging SE")

#----------------------#
#kriging with gstat
#----------------------#
new.grid.1m.gstat <- expand.grid(x=seq(0,max(matrix$x),.5), y=seq(0,max(matrix$y),.5))#need labels for coords

#convert to sp object
gridded(new.grid.1m.gstat) = ~x + y

krig.gstat <- krige(Height ~ 1, gstat.veg, new.grid.1m.gstat, model = sph.gstat)

#Plot
#image(krig.gstat, main="krigged estimates-gstat")
#image(krig.geoR.exp, main="krigged estimates-geoR")

cor(krig.geoR.exp$predict, krig.gstat$var1.pred)

#Inverse distance weighting in gstat
idw.gstat <- idw(Height ~ 1, gstat.veg, new.grid.1m.gstat)#idp-the power used; default power=2


#Correlation between idw, kriging in gstat, kriging in geoR
round(cor(cbind(geoR.exp=krig.geoR.exp$predict,
                gstat.exp=krig.gstat$var1.pred,
                gstat.idw=idw.gstat$var1.pred)), 3)

################################################################################
#code for question 9

#download elevation code
#change the directory as needed
srtm <- raster("~/Desktop/GitHub/ECOL_620/data/data_for_lab5/srtm_15_05.tif")
e <- extent(-107.5, -107, 37.5, 38)
srtm=crop(srtm, e)
plot(srtm)
#aggregate to reduce computation
srtm = aggregate(srtm, fact=10, fun=mean) #mean value
plot(srtm)

#extract 500 points of elevation, randomly drawn
set.seed(2020)
srtm_pts=as.data.frame(sampleRandom(srtm, 500, xy=T)) #xy = T means keep x and y coords
colnames(srtm_pts)[3]="elevation"
head(srtm_pts)
plot(srtm)
points(srtm_pts$x,srtm_pts$y)

#Create grid with intervals of 0.008333333 degree # prdict kriging back onto this 
extent(srtm)
res(srtm)
e <- extent(srtm)
new.grid<- raster(e, resolution=res(srtm))
new.grid=as.data.frame(new.grid, xy=TRUE)
new.grid=new.grid[,1:2]
